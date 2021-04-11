using Parameters

using Pkg
# Pkg.add("Statistics")
using Interpolations

# For testing
using Random
using Plots
using Distributions
using Statistics

##

# Structure to hold all parameters and the setup functions of the problem

@with_kw mutable struct ModelDiscrete

    # Discount factor
    beta0 = 0.9

    # Relative Risk Aversion
    gamma = 1.0

    # Use CRRA?
    crra = 1 # 0 means Epstein-Zin

    r = 0.01
    R = 1+r

    # Build grids
    na = 50
    a_bc = 0
    a_max = 200
    agrid::Array{Float64,1} = LinRange(a_bc,a_max,na).^(0.5) # Needs to be properly built

    # Used as an indicator during EGP (stored here for inplace updating)
    borrow_constr = zeros(na)

    # Stores interpolator used in EGP
    tmp_itp = NaN

    nyF = 2
    yFgrid::Array{Float64,1} = LinRange(1.0,1.0,nyF).^(0.5)
    nyP = 2
    yPgrid::Array{Float64,1} = LinRange(0,0.1,nyP).^(0.5)
    nyT = 20
    yTgrid::Array{Float64,1} = LinRange(0,0.1,nyT).^(0.5)

    nz = 2
    zgrid = LinRange(0,nz,nz)

    # Build augmented grids (Fix this!)
    aagrid = repeat(reshape(agrid, (na,1,1,1,1)), 1, nyF, nyP, nyT, nz)

    yyFgrid = repeat(reshape(yFgrid, (1,nyF,1,1,1)), na, 1, nyP, nyT, nz)

    yyPgrid = repeat(reshape(yPgrid, (1,1,nyP,1,1)), na, nyF, 1, nyT, nz)

    yyTgrid = repeat(reshape(yTgrid, (1,1,1,nyT,1)), na, nyF, nyP, 1, nz)

    zzgrid = repeat(reshape(zgrid, (1,1,1,1,nz)), na, nyF, nyP, nyT, 1)

    # income transitions (needs to be set up)
    income_trans = [zeros(na,nyF,nyP,nyT,nz)/(na*nyF*nyP*nyT*nz) for a in 1:na, yF in 1:nyF, yP in 1:nyP, yT in 1:nyT, z in 1:nz]

    # Consumption as function of assets and income today
    con = zeros(na,nyF,nyP,nyT,nz)

    # Consumption as function of assets tomorrow and income today
    con_euler = zeros(na,nyF,nyP,nyT,nz)

    # Asset choice as function of assets tomorrow and income today
    a_euler = zeros(na,nyF,nyP,nyT,nz)

    # choice of assets held (saved)
    a_tom = zeros(na,nyF,nyP,nyT,nz)

    # u'(c)
    up = zeros(na,nyF,nyP,nyT,nz)

    # [u']^{-1}(c)
    upinv = zeros(na,nyF,nyP,nyT,nz) # Not needed anymore

    # MUC
    muc = zeros(na,nyF,nyP,nyT,nz)

    # Expected MUC
    emuc = zeros(na,nyF,nyP,nyT,nz)

    egp_maxiter = 1000
    egp_tol = 1e-6

end

# Marginal utility
function up(m)
    # println("up called")
    if m.crra == 1
        m.up = m.con.^(-m.gamma)
    else
        # Do Epstein-Zin
        m.up = NaN
    end
end


# Expected MUC
function emuc(m)
    # Should add threading
    for a in 1:m.na
        for yF in 1:m.nyF
            for yP in 1:m.nyP
                for yT in 1:m.nyT
                    for z in 1:m.nz
                        m.emuc[a,yF,yP,yT,z] = sum(m.up .* m.income_trans[a,yF,yP,yT,z])
                    end
                end
            end
        end
    end
end

function muc(m)
    m.muc = m.beta0 .* m.R .* m.emuc
end

# function upinv(p)
#     if m.crra == 1
#         m.upinv = m.muc.^(-1/m.gamma)
#     else
#         m.upinv = NaN
#     end
# end

function update_con_euler(m)
    if m.crra == 1
        m.con_euler = m.muc.^(-1/m.gamma)
    else
        m.con_euler = NaN
    end
end

function update_a_euler(m)
    m.a_euler = (m.con_euler + m.aagrid - m.yyFgrid - m.yyPgrid - m.yyTgrid)/m.R
end

function interp_a_tom(m)
    for yF in 1:m.nyF
        for yP in 1:m.nyP
            for yT in 1:m.nyT
                for z in 1:m.nz

                    m.borrow_constr = (m.agrid .< m.a_euler[1,yF,yP,yT,z])

                    m.tmp_itp = interpolate((m.a_euler[:,yF,yP,yT,z],), m.agrid, Gridded(Linear()))

                    m.a_tom[:,yF,yP,yT,z] = (1 .- m.borrow_constr) .* extrapolate(m.tmp_itp, Line()).(m.agrid) .+ m.borrow_constr .* m.agrid[1]

                    # m.a_tom[:,yF,yP,yT,z] = (m.a_tom[:,yF,yP,yT,z] .>= m.a_bc) .* m.a_tom[:,yF,yP,yT,z] .+ (m.a_tom[:,yF,yP,yT,z] .< m.a_bc) .* m.agrid[1]

                    # m.a_tom[:,yF,yP,yT,z] = m.borrow_constr .* extrapolate(m.tmp_itp, Line()).(m.agrid) .+ (1 .- m.borrow_constr) .* m.agrid[1]
                end
            end
        end
    end
end

function update_con(m)
    m.con = m.R*m.aagrid + m.yyFgrid + m.yyPgrid + m.yyTgrid - m.a_tom
end


function setup_income(m)
    m.income_trans *= 0.0
    # Add threading?
    for a in 1:m.na
        for yF in 1:m.nyF
            for yP in 1:m.nyP
                for yT in 1:m.nyT
                    for z in 1:m.nz
                        for a2 in 1:m.na
                            if a == a2
                                # println("size(ones(yF,yP,yT,z)./sum(ones(yF,yP,yT,z))): ", size(ones(m.nyF,m.nyP,m.nyT,m.nz)./sum(ones(m.nyF,m.nyP,m.nyT,m.nz))))
                                m.income_trans[a,yF,yP,yT,z][a2,:,:,:,:] = ones(m.nyF,m.nyP,m.nyT,m.nz)./sum(ones(m.nyF,m.nyP,m.nyT,m.nz))
                            end
                        end
                    end
                end
            end
        end
    end
end



##

# p0 = Params()
# println(p0)

# @unpack na,nyF,nyP,nyT,nz = p0
# c0 = exm.(randn(na,nyF,nyP,nyT,nz))

# println("1. maximum(p0.up): ", maximum(p0.up), ", size(p0.up): ", size(p0.up))
# up(p0, c0[1,1,1,1,1])
# println("2. maximum(p0.up): ", maximum(p0.up), ", size(p0.up): ", size(p0.up))
# up(p0, c0)
# println("3. maximum(p0.up): ", maximum(p0.up), ", size(p0.up): ", size(p0.up))
# emuc(p0,c0)
# p0.emuc
# p0.up = 5
# println("2. maximum(p0.up): ", maximum(p0.up))
# println("1. r: ", p0.r, ", R: ", p0.R)
# p0.r = 0.5
# println("2. r: ", p0.r, ", R: ", p0.R)

# println(p0.agrid)

# itp = interpolate(([1, 2, 3],), [4, 5, 6], Gridded(Linear()))
# # itp = interpolate([1 2 3; 4 5 6], BSpline(Linear()))
# itm.([2 2.3 2.4 2.7 3])
# extp = extrapolate(itp, Line())
# extm.([0.5 1 1.5 2 3 3.5])