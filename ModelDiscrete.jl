using Parameters

using Pkg
# Pkg.add("DelimitedFiles")
using Interpolations

# For testing
using Random
using Plots
using Distributions
using Statistics
using DelimitedFiles

using LinearAlgebra

using Roots

##

# Structure to hold all parameters and the setup functions of the problem

@with_kw mutable struct ModelDiscrete

    # Name of ModelDiscrete
    name = "default"

    # Description
    descr = "default model"

    # Mean annual income in dollars
    numeraire_in_dollars = 72000

    # Discount factor (initial)
    beta0 = 0.95

    # Discount factor (bounds for search)

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
    yFgrid::Array{Float64,1} = exp.(LinRange(0.0,.0,nyF))
    yFtrans = Matrix{Float64}(I,nyF,nyF)

    # Persistent
    nyP = 4
    yP_sd = 0.1
    yP_rho = 0.0
    yPgrid::Array{Float64,1} = exp.(rouwenhorst(nyP, -(yP_sd^2)/2,yP_sd, yP_rho)[1])
    yPtrans = rouwenhorst(nyP, -(yP_sd^2)/2,yP_sd, yP_rho)[2]

    # Transitory
    nyT = 4
    yT_sd = 0.2
    yTgrid::Array{Float64,1} = exp.(match_normal(nyT, 0, yT_sd)[2])
    yTtrans = repeat(reshape(match_normal(nyT, 0, yT_sd)[3], (1,nyT)), nyT)

    # Discount factor heterogeneity
    nbh = 2
    bh_sd = 0
    bhgrid = LinRange(0,nbh,nbh)
    bhtrans = repeat(reshape(match_normal(nbh, 0, bh_sd)[3], (1,nbh)), nbh)

    # Build augmented grids (Fix this!)
    aagrid = repeat(reshape(agrid, (na,1,1,1,1)), 1, nyF, nyP, nyT, nbh)

    yyFgrid = repeat(reshape(yFgrid, (1,nyF,1,1,1)), na, 1, nyP, nyT, nbh)

    yyPgrid = repeat(reshape(yPgrid, (1,1,nyP,1,1)), na, nyF, 1, nyT, nbh)

    yyTgrid = repeat(reshape(yTgrid, (1,1,1,nyT,1)), na, nyF, nyP, 1, nbh)

    bbhgrid = repeat(reshape(bhgrid, (1,1,1,1,nbh)), na, nyF, nyP, nyT, 1)

    # Income transitions (needs to be set up)
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
    target_mean_wealth = 10.2

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
                        # m.income_trans[a,yF,yP,yT,bh][a,:,:,:,:] = ones(m.nyF,m.nyP,m.nyT,m.nbh)./sum(ones(m.nyF,m.nyP,m.nyT,m.nbh))
                        
                        # Fixed
                        m.income_trans[a,yF,yP,yT,bh][a,:,:,:,:] = repeat(reshape(m.yFtrans[yF,:], (m.nyF,1,1,1)), 1, m.nyP, m.nyT, m.nbh)
                        # Persistent
                        m.income_trans[a,yF,yP,yT,bh][a,:,:,:,:] .*= repeat(reshape(m.yPtrans[yP,:], (1,m.nyP,1,1)), m.nyF, 1, m.nyT, m.nbh)
                        # Transitory
                        m.income_trans[a,yF,yP,yT,bh][a,:,:,:,:] .*= repeat(reshape(m.yTtrans[yT,:], (1,1,m.nyT,1)), m.nyF, m.nyP, 1, m.nbh)

                        # beta
                        m.income_trans[a,yF,yP,yT,bh][a,:,:,:,:] .*= repeat(reshape(m.bhtrans[bh,:], (1,1,1,m.nbh)), m.nyF, m.nyP, m.nyT, 1)
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

# Creates equispaced approximation to normal distribution
# n - number of points
# mu - mean
# sigma - standard deviation
# width - multiple of standard deviation for width of grid
# returns:
    # f - error in approximation
    # x - location of points
    # p - probabilities
function discrete_normal(n, mu, sigma, width)
    if n > 2
        x = collect(LinRange(mu - width*sigma, mu + width*sigma, n))
        # discretize CDF trapezoid-wise
        p = zeros(n)
        p[1] = cdf(Normal(mu,sigma), x[1] + (x[2] - x[1])/2)
        for i in 2:n-1
            p[i] = cdf(Normal(mu,sigma), x[i] + (x[i+1] - x[i])/2) - cdf(Normal(mu,sigma), x[i] - (x[i] - x[i-1])/2)
        end
        p[n] = 1 - sum(p[1:n-1])
    elseif n == 2
        x = collect(LinRange(mu - width*sigma, mu + width*sigma, n))
        # split the mass
        p = 0.5 * ones(n)
    else
        x = [mu]
        p = [1.0]
    end
    
    Ex = x'*p
    SDx = sqrt.((x'.^2)*p .- Ex^2)
    
    # Error between realized sd and desired
    f = SDx .- sigma
    return [f, x, p]
end


# Probably move these to their own .jl files?
# Returns normal approximation with optimal width
function match_normal(n, mu, sigma, width_guess=0.5)
    width = find_zero(x -> discrete_normal(n,mu,sigma,x),width_guess)
    temp, grid, dist = discrete_normal(n, mu, sigma, width)
end


# Rouwenhorst discretization of AR(1)
function rouwenhorst(n, mu, sigma, rho)
    # Get width of grid
    width = sqrt((n-1)* sigma^2/(1-rho^2)) # This is correct, though both Greg and Kopecky have typos in their note
    
    # Equi-space grid
    grid = LinRange(mu - width, mu + width, n)

    # Initialize for 2-case
    p = (1 + rho)/2
    trans = [p 1-p; 1-p p]

    # Recursively find n-case
    if n > 2
        for i in 2:n-1
            # Per Rouwenhorst
            trans = p.*[trans zeros(i); zeros(i)' 0] + (1-p).*[zeros(i) trans; 0 zeros(i)'] + (1-p).*[zeros(i)' 0; trans zeros(i)] + p .* [0 zeros(i)'; zeros(i) trans]
            # Row sums should be 2 for {2,...,n-1}, and 1 for {1,n}
            trans ./= sum(trans, dims=2)
        end
    end

    return grid, trans
end

# row_test = rouwenhorst(2,2,0.1,0.8)
# println("row grid: ", row_test[1])
# println("row trans: ", row_test[2])
# row_test[2]
# sum(row_test[2], dims=2)

# readdlm("quarterly_b_yPtrans.txt")

# m0 = ModelDiscrete()

# setup_income(m0)

# size(m0.income_trans)

# m0.yFgrid
# m0.yFtrans

# m0.yPgrid
# m0.yPtrans




##

# match_normal(5, 2, 0)

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