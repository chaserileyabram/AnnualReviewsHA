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

    # Discount factor (initial)
    beta0 = 0.9

    # Relative Risk Aversion
    gamma = 1.0

    # Use CRRA?
    crra = 1 # 0 means Epstein-Zin

    # Interest Rate
    r = 0.01
    R = 1+r

    # Build grids
    na = 10
    # Borrowing constraint
    a_bc = 0
    # Max assets
    a_max = 200
    # Shape for power spacing
    a_shape = 0.15
    # Asset grid (not set up)
    agrid::Array{Float64,1} = LinRange(0,1,na)

    # Indicate whether asset grid has been set up
    grid_setup = false

    # Used as an indicator during EGP (stored here for inplace updating)
    borrow_constr =  zeros(na)

    # Stores interpolator used in EGP
    tmp_itp = NaN

    # Income
    # Fixed
    nyF = 2
    yFgrid::Array{Float64,1} = LinRange(1.0,1.0,nyF)
    # Persisten
    nyP = 2
    yPgrid::Array{Float64,1} = LinRange(0.9,1.1,nyP)
    # Transitory
    nyT = 2
    yTgrid::Array{Float64,1} = LinRange(0.5,1.5,nyT)

    # Discount factor heterogeneity
    nbh = 2
    bhgrid = LinRange(0,nbh,nbh)

    # Build augmented grids (Fix this!)
    aagrid = repeat(reshape(agrid, (na,1,1,1,1)), 1, nyF, nyP, nyT, nbh)

    yyFgrid = repeat(reshape(yFgrid, (1,nyF,1,1,1)), na, 1, nyP, nyT, nbh)

    yyPgrid = repeat(reshape(yPgrid, (1,1,nyP,1,1)), na, nyF, 1, nyT, nbh)

    yyTgrid = repeat(reshape(yTgrid, (1,1,1,nyT,1)), na, nyF, nyP, 1, nbh)

    bbhgrid = repeat(reshape(bhgrid, (1,1,1,1,nbh)), na, nyF, nyP, nyT, 1)

    # income transitions (needs to be set up)
    # This is a' as a function of a', so no movement in a dim
    income_trans = [zeros(na,nyF,nyP,nyT,nbh)/(na*nyF*nyP*nyT*nbh) for a in 1:na, yF in 1:nyF, yP in 1:nyP, yT in 1:nyT, bh in 1:nbh]

    # Same shape as income, but now account for policy in a (and/or taxes)
    # a' as a function of a
    inter_trans = [zeros(na,nyF,nyP,nyT,nbh)/(na*nyF*nyP*nyT*nbh) for a in 1:na, yF in 1:nyF, yP in 1:nyP, yT in 1:nyT, bh in 1:nbh]

    # Used for interpolating in building inter_trans
    # lower bound
    lb = 0.0
    # upper bound
    ub = 1.0
    # mixing weight on lower bound
    lbubmix = 0.5

    # Full transition (requires policy solved)
    full_trans_dim = na*nyF*nyP*nyT*nbh
    full_trans =  zeros(full_trans_dim, full_trans_dim)

    # Consumption as function of assets and income today
    con =  zeros(na,nyF,nyP,nyT,nbh)

    # Consumption as function of assets tomorrow and income today
    con_euler =  zeros(na,nyF,nyP,nyT,nbh)

    # Asset choice as function of assets tomorrow and income today
    a_euler =  zeros(na,nyF,nyP,nyT,nbh)

    # choice of assets held (saved)
    a_tom =  zeros(na,nyF,nyP,nyT,nbh)

    # u'(c)
    up =  zeros(na,nyF,nyP,nyT,nbh)

    # [u']^{-1}(c)
    upinv =  zeros(na,nyF,nyP,nyT,nbh) # Not needed anymore

    # MUC
    muc =  zeros(na,nyF,nyP,nyT,nbh)

    # Expected MUC
    emuc =  zeros(na,nyF,nyP,nyT,nbh)

    # Stationary Distribution
    statdist = ones(na,nyF,nyP,nyT,nbh)./full_trans_dim
    statdistvec = ones(full_trans_dim)./full_trans_dim

    # Distribution over just assets
    adist = ones(na)/na

    # Iteration pars on EGP
    egp_maxiter = 1000
    egp_tol = 1e-6

    # Iteration pars on stationary distribution
    statdist_maxiter = 1000
    statdist_tol = 1e-9

    # Targets
    target_mean_wealth = 0.9

end

# Marginal utility
function up(m)
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
                    for bh in 1:m.nbh
                        m.emuc[a,yF,yP,yT,bh] = sum(m.up .* m.income_trans[a,yF,yP,yT,bh])
                    end
                end
            end
        end
    end
end

# Marginal utility of consumption
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

# Update consumption used for Euler
function update_con_euler(m)
    if m.crra == 1
        m.con_euler = m.muc.^(-1/m.gamma)
    else
        m.con_euler = NaN
    end
end

# Update assets used for Euler
function update_a_euler(m)
    m.a_euler = (m.con_euler + m.aagrid - m.yyFgrid .* m.yyPgrid .* m.yyTgrid)/m.R
end

# Use interpolation to "flip" savings function
# a(a') -> a'(a)
function interp_a_tom(m)
    for yF in 1:m.nyF
        for yP in 1:m.nyP
            for yT in 1:m.nyT
                for bh in 1:m.nbh
                    # Find where borrowing constrained
                    m.borrow_constr = (m.agrid .< m.a_euler[1,yF,yP,yT,bh])
                    # Build interpolating function (just linear)
                    m.tmp_itp = interpolate((m.a_euler[:,yF,yP,yT,bh],), m.agrid, Gridded(Linear()))
                    # Use interpolation to invert
                    m.a_tom[:,yF,yP,yT,bh] = (1 .- m.borrow_constr) .* extrapolate(m.tmp_itp, Line()).(m.agrid) .+ m.borrow_constr .* m.agrid[1]
                end
            end
        end
    end
end

# Update consumption (need a_tom updated)
function update_con(m)
    m.con = m.R*m.aagrid + m.yyFgrid .* m.yyPgrid .* m.yyTgrid - m.a_tom
end

# Set up income process
# currently iid and dumb
function setup_income(m)
    # Zero out
    m.income_trans *= 0.0
    # Add threading?
    for a in 1:m.na
        for yF in 1:m.nyF
            for yP in 1:m.nyP
                for yT in 1:m.nyT
                    for bh in 1:m.nbh
                        m.income_trans[a,yF,yP,yT,bh][a,:,:,:,:] = ones(m.nyF,m.nyP,m.nyT,m.nbh)./sum(ones(m.nyF,m.nyP,m.nyT,m.nbh))
                    end
                end
            end
        end
    end
end

# Setup grids with power spacing
function setup_power_grids(m)
    # Start with [0,1]
    m.agrid = LinRange(0,1,m.na)
    # Powerspace
    m.agrid = m.agrid.^(1/m.a_shape)
    # Scale to bounds
    m.agrid = m.a_bc .+ (m.a_max - m.a_bc) .* m.agrid

    # Apply to aagrid
    m.aagrid = repeat(reshape(m.agrid, (m.na,1,1,1,1)), 1, m.nyF, m.nyP, m.nyT, m.nbh)

    # Indicate that grids have been setup
    m.grid_setup = true
end


##

# p0 = Params()
# println(p0)

# @unpack na,nyF,nyP,nyT,nbh = p0
# c0 = exm.(randn(na,nyF,nyP,nyT,nbh))

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