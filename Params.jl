using Parameters
using Random # only for testing

using Pkg
# Pkg.add("Interpolations")
using Interpolations

##

# Structure to hold all parameters and the setup functions of the problem

@with_kw mutable struct Params

    # Discount factor
    beta0 = 0.95

    # Relative Risk Aversion
    gamma = 1.0

    # Use CRRA?
    crra = 1 # 0 means Epstein-Zin

    r = 0.01
    R = 1+r

    # Build grids
    na = 5
    a_bc = 0
    a_max = 200
    agrid = LinRange(0,na,na) # Needs to be properly built

    nyF = 2
    yFgrid = LinRange(0,nyF,nyF)
    nyP = 3
    yPgrid = LinRange(0,nyP,nyP)
    nyT = 4
    yTgrid = LinRange(0,nyT,nyT)

    nz = 3
    zgrid = LinRange(0,nz,nz)

    # Build augmented grids (Fix this!)
    aagrid = repeat(reshape(agrid, (na,1,1,1,1)), 1, nyF, nyP, nyT, nz)

    yyFgrid = repeat(reshape(yFgrid, (1,nyF,1,1,1)), na, 1, nyP, nyT, nz)

    yyPgrid = repeat(reshape(yPgrid, (1,1,nyP,1,1)), na, nyF, 1, nyT, nz)

    yyTgrid = repeat(reshape(yTgrid, (1,1,1,nyT,1)), na, nyF, nyP, 1, nz)

    zzgrid = repeat(reshape(zgrid, (1,1,1,1,nz)), na, nyF, nyP, nyT, 1)



    # income transitions (needs to be set up)
    income_trans = [ones(na,nyF,nyP,nyT,nz)/(na*nyF*nyP*nyT*nz) for a in 1:na, yF in 1:nyF, yP in 1:nyP, yT in 1:nyT, z in 1:nz]

    # Consumption as function of assets tomorrow and income today
    con_euler = zeros(na,nyF,nyP,nyT,nz)

    # Asset choice as function of assets tomorrow and income today
    a_euler = zeros(na,nyF,nyP,nyT,nz)

    # choice of assets held (saved)
    a_tom = zeros(na,nyF,nyP,nyT,nz)

    # u'(c)
    up = zeros(na,nyF,nyP,nyT,nz)

    # [u']^{-1}(c)
    upinv = zeros(na,nyF,nyP,nyT,nz)

    # MUC
    muc = zeros(na,nyF,nyP,nyT,nz)

    # Expected MUC
    emuc = zeros(na,nyF,nyP,nyT,nz)

    egp_maxiter = 100
    egp_tol = 1e-6

end

# Marginal utility
function up(p)
    if p.crra == 1
        p.up = p.con.^(-p.gamma)
    else
        # Do Epstein-Zin
        p.up = NaN
    end
end

function upinv(p)
    if p.crra == 1
        p.upinv = p.muc.^(-1/p.gamma)
    else
        p.upinv = NaN
    end
end


# Expected MUC
function emuc(p)
    # Should add threading
    for a in 1:p.na
        for yF in 1:p.nyF
            for yP in 1:p.nyP
                for yT in 1:nyT
                    for z in 1:p.nz
                        p.emuc[a,yF,yP,yT,z] = sum(p.up .* p.income_trans[a,yF,yP,yT,z])
                    end
                end
            end
        end
    end
end

function muc(p)
    p.muc = p.beta0*p.R*p.emuc
end

function update_con_euler(p)
    p.con_euler = upinv(p,p.muc)
end

function update_a_euler(p)
    p.a_euler = 1/p.R*(p.con + p.aagrid - p.yyFgrid - p.yyPgrid - p.yyTgrid)
end

function interp_a_tom(p)
    
    for yF in 1:p.nyF
        for yP in 1:p.nyP
            for yT in 1:nyT
                for z in 1:p.nz
                    borrow_constr = (p.agrid .< p.a_euler[1,yF,yP,yT,z])
                    p.a_tom[:,yF,yP,yT,z] = borrow_constr .* interpolate(p.a_euler[:,yF,yP,yT,z], p.agrid).(p.agrid) .+ (1 .- borrow_constr) .* agrid[1]
                end
            end
        end
    end
end

function update_con(p)
    p.con = p.R*p.aagrid + p.yyFgrid + p.yyPgrid + p.yyTgrid - p.a_tom
end








##

# p0 = Params()
# println(p0)

# @unpack na,nyF,nyP,nyT,nz = p0
# c0 = exp.(randn(na,nyF,nyP,nyT,nz))

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
