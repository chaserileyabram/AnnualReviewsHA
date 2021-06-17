# Struct to hold model parameters for 2-asset continuous model

# Chase Abram

using Parameters
using LinearAlgebra
using Roots

using NLsolve
using ForwardDiff

using SparseArrays

using Distributions
##

@with_kw mutable struct ModelContinuous

    # Name of ModelDiscrete
    name = "default"

    # Description
    descr = "default continuous model"

    # Mean annual income in dollars
    numeraire_in_dollars = 72000

    # Discount rate (initial)
    rho = 0.05/4

    # Discount factor (bounds for search)

    # Relative Risk Aversion
    gamma = 1.0

    # Use CRRA?
    crra = 1.0 # 0 means stochastic differential utility (not implemented)

    # Interest Rates
    # Liquid
    rb = 0.01/4
    Rb = 1+rb

    # Liquid
    ra = 0.06/4
    Ra = 1+ra

    # Death rate
    zeta = 0.0

    # Build grids

    # Liquid
    nb = 30
    # Borrowing constraint
    b_bc = 0
    # Max assets
    b_max = 20
    # Shape for power spacing
    b_shape = 0.15
    # Asset grid (not set up)
    bgrid::Array{Float64,1} = power_grid(nb,b_shape,b_bc,b_max)

    # Illiquid
    na = 31
    # Borrowing constraint
    a_bc = 0
    # Max assets
    a_max = 500
    # Shape for power spacing
    a_shape = 0.15
    # Asset grid (not set up)
    agrid::Array{Float64,1} = power_grid(na,a_shape,a_bc,a_max)

    # Flow adjustment cost parameters
    
    # chi(d) = chi_0 *|d| + chi_1/(1 + chi_2) *|d/a|^(1 + chi_2)
    chi_0 = 0.0
    chi_1 = 5.0
    chi_2 = 1.0

    a_lb = 0.25
    dep_max = 100
    chid_max = 0.0
    a_stable = 0.0

    # Indicate whether asset grid has been set up
    # Should set up all grids in here?
    # grid_setup = false

    # Used as an indicator during EGP (stored here for inplace updating)
    # borrow_constr =  zeros(na)

    # Stores interpolator used in EGP
    # tmp_itp = NaN

    # Income
    income_arrival = 0.25

    # Fixed
    nyF = 2
    yFgrid::Array{Float64,1} = exp.(LinRange(0.0,0.0,nyF))
    yFtrans = Matrix{Float64}(I,nyF,nyF)

    # Persistent
    nyP = 3
    yP_sd = 0 #sqrt(0.01080)
    yP_rho = 0.9881
    yPgrid::Array{Float64,1} = exp.(rouwenhorst(nyP, -(yP_sd^2)/2,yP_sd, yP_rho)[1])
    yPtrans = rouwenhorst(nyP, -(yP_sd^2)/2,yP_sd, yP_rho)[2]

    # Transitory
    nyT = 2
    yT_sd = 0 #sqrt(0.2087)
    yTgrid::Array{Float64,1} = exp.(match_normal(nyT, 0, yT_sd)[2])
    yTtrans = repeat(reshape(match_normal(nyT, 0, yT_sd)[3], (1,nyT)), nyT)

    # Discount factor heterogeneity
    nrh = 2
    rh_sd = 0
    rhgrid = LinRange(0,nrh,nrh)
    rhtrans = repeat(reshape(match_normal(nrh, 0, rh_sd)[3], (1,nrh)), nrh)

    # Build augmented grids (Fix this for more het, if needed)
    bbgrid = repeat(reshape(bgrid, (nb,1,1,1,1,1)), 1, na, nyF, nyP, nyT, nrh)

    aagrid = repeat(reshape(agrid, (1,na,1,1,1,1)), nb, 1, nyF, nyP, nyT, nrh)

    yyFgrid = repeat(reshape(yFgrid, (1,1,nyF,1,1,1)), nb, na, 1, nyP, nyT, nrh)

    yyPgrid = repeat(reshape(yPgrid, (1,1,1,nyP,1,1)), nb, na, nyF, 1, nyT, nrh)

    yyTgrid = repeat(reshape(yTgrid, (1,1,1,1,nyT,1)), nb, na, nyF, nyP, 1, nrh)

    rrhgrid = repeat(reshape(rhgrid, (1,1,1,1,1,nrh)), nb, na, nyF, nyP, nyT, 1)


    # Some of this was for EGP
    # Income transitions (needs to be set up)
    # # This is a' as a function of a', so no movement in a dim
    # income_trans = [zeros(nb,na,nyF,nyP,nyT,nrh)/(nb*na*nyF*nyP*nyT*nrh) for b in 1:nb, a in 1:na, yF in 1:nyF, yP in 1:nyP, yT in 1:nyT, rh in 1:nrh]

    # # Same shape as income, but now account for policy in a (and/or taxes)
    # # a' as a function of a
    # inter_trans = [zeros(na,nyF,nyP,nyT,nrh)/(na*nyF*nyP*nyT*nrh) for a in 1:na, yF in 1:nyF, yP in 1:nyP, yT in 1:nyT, rh in 1:nrh]

    # # Used for interpolating in building inter_trans
    # # lower bound
    # lb = 0.0
    # # upper bound
    # ub = 1.0
    # # mixing weight on lower bound
    # lbubmix = 0.5

    # # Full transition (requires policy solved)
    full_trans_dim = nb*na*nyF*nyP*nyT*nrh
    # full_trans =  zeros(full_trans_dim, full_trans_dim)

    # Value function
    V = zeros(nb,na,nyF,nyP,nyT,nrh)

    # Finite differences
    dVFb = fill(NaN,nb,na,nyF,nyP,nyT,nrh)
    dVBb = fill(NaN,nb,na,nyF,nyP,nyT,nrh)
    dVbmin = 1e-8

    dVFa = fill(NaN,nb,na,nyF,nyP,nyT,nrh)
    dVBa = fill(NaN,nb,na,nyF,nyP,nyT,nrh)
    dVamin = 0.0


    Vvec = zeros(nb*na*nyF*nyP*nyT*nrh,1)
    Vvecnew = zeros(nb*na*nyF*nyP*nyT*nrh,1)

    # Consumption as function of assets and income today
    conF =  zeros(nb,na,nyF,nyP,nyT,nrh)
    conB =  zeros(nb,na,nyF,nyP,nyT,nrh)
    con0 =  zeros(nb,na,nyF,nyP,nyT,nrh)
    con =  zeros(nb,na,nyF,nyP,nyT,nrh)

    # Savings (consumption part)
    savcF =  zeros(nb,na,nyF,nyP,nyT,nrh)
    savcB =  zeros(nb,na,nyF,nyP,nyT,nrh)
    savc0 =  zeros(nb,na,nyF,nyP,nyT,nrh)
    savc =  zeros(nb,na,nyF,nyP,nyT,nrh)

    # Savings (deposit part)
    # savdF =  zeros(nb,na,nyF,nyP,nyT,nrh)
    # savdB =  zeros(nb,na,nyF,nyP,nyT,nrh)
    # savd =  zeros(nb,na,nyF,nyP,nyT,nrh)

    # Hamiltonians
    HcF = zeros(nb,na,nyF,nyP,nyT,nrh)
    HcB = zeros(nb,na,nyF,nyP,nyT,nrh)
    Hc0 = zeros(nb,na,nyF,nyP,nyT,nrh)

    # Upwind direction indicators
    IcF = zeros(nb,na,nyF,nyP,nyT,nrh)
    IcB = zeros(nb,na,nyF,nyP,nyT,nrh)
    Ic0 = zeros(nb,na,nyF,nyP,nyT,nrh)


    # Consumption as function of assets tomorrow and income today
    # con_euler =  zeros(na,nyF,nyP,nyT,nrh)

    # Asset choice as function of assets tomorrow and income today
    # a_euler =  zeros(na,nyF,nyP,nyT,nrh)

    # choice of assets held (saved)
    # a_tom =  zeros(na,nyF,nyP,nyT,nrh)

    # Deposits
    depFbFa = zeros(nb,na,nyF,nyP,nyT,nrh)
    depFbBa = zeros(nb,na,nyF,nyP,nyT,nrh)
    depBbFa = zeros(nb,na,nyF,nyP,nyT,nrh)
    depBbBa = zeros(nb,na,nyF,nyP,nyT,nrh)

    # depF = zeros(nb,na,nyF,nyP,nyT,nrh)
    # depB = zeros(nb,na,nyF,nyP,nyT,nrh)
    dep = zeros(nb,na,nyF,nyP,nyT,nrh)

    savdFF =  zeros(nb,na,nyF,nyP,nyT,nrh)
    savdFB =  zeros(nb,na,nyF,nyP,nyT,nrh)
    savdBF =  zeros(nb,na,nyF,nyP,nyT,nrh)
    savdBB =  zeros(nb,na,nyF,nyP,nyT,nrh)
    savd =  zeros(nb,na,nyF,nyP,nyT,nrh)

    # Hamiltonians
    HdFF = zeros(nb,na,nyF,nyP,nyT,nrh)
    HdFB = zeros(nb,na,nyF,nyP,nyT,nrh)
    HdBF = zeros(nb,na,nyF,nyP,nyT,nrh)
    HdBB = zeros(nb,na,nyF,nyP,nyT,nrh)

    # Upwind direction indicators
    IdFF = zeros(nb,na,nyF,nyP,nyT,nrh)
    IdFB = zeros(nb,na,nyF,nyP,nyT,nrh)
    IdBF = zeros(nb,na,nyF,nyP,nyT,nrh)
    IdBB = zeros(nb,na,nyF,nyP,nyT,nrh)


    # Choices  (no income, death, etc.)
    A = spzeros(full_trans_dim, full_trans_dim)

    # Income transitions
    AL = spzeros(full_trans_dim, full_trans_dim)

    # Death and Birth transitions
    AD = spzeros(full_trans_dim, full_trans_dim)

    # u'(c)
    # up =  zeros(na,nyF,nyP,nyT,nrh)

    # [u']^{-1}(c)
    # upinv =  zeros(na,nyF,nyP,nyT,nrh) # Not needed anymore

    # MUC
    # muc =  zeros(na,nyF,nyP,nyT,nrh)

    # Expected MUC
    # emuc =  zeros(na,nyF,nyP,nyT,nrh)

    # Stationary Distribution
    statdist = ones(na,nyF,nyP,nyT,nrh)./full_trans_dim
    statdistvec = ones(full_trans_dim)./full_trans_dim

    # Distribution over just liquid assets
    bdist = ones(na)/na

    # Distribution over just illiquid assets
    adist = ones(na)/na

    # Iteration pars on EGP
    hjb_maxiter = 1000
    hjb_tol = 1e-9
    hjb_delta = 1e6

    # Iteration pars on stationary distribution
    kfe_maxiter = 1000
    kfe_tol = 1e-9

    # Targets
    target_mean_wealth = 4.1
    target_mean_liquid = 0.562

end

# Functions to setup should go here, 
# so that creating the model sets it up to be solved immediately

function uc(m,z)
    if m.gamma == 1
        return log.(z)
    else
        return (z.^(1-m.gamma) - 1)./(1-m.gamma)
    end
end

function ucp_inv(m,z)
    return z.^(-1/m.gamma)
end

# chi(d) = chi_0 *|d| + chi_1/(1 + chi_2) *|d/a|^(1 + chi_2)
function chi(m,d,a)
    return m.chi_0*abs.(d) + m.chi_1/(1 + m.chi_2).*abs.(d).^(1 + m.chi_2) ./ (max.(a,m.a_lb).^m.chi_2)
end

function chi_d(m,d,a)
    # if d > 0
    #     return m.chi_0 .+ m.chi_1.*abs.(d./max.(a,m.a_lb)).^m.chi_2
    # else
    #     return -m.chi_0 .- m.chi_1.*abs.(d./max.(a,m.a_lb)).^m.chi_2
    # end

    return sign(d).*(m.chi_0 .+ m.chi_1.*abs.(d./max.(a,m.a_lb)).^m.chi_2)
end

function chi_inv(m,chid,a)

    m.a_stable = max(a, m.a_lb)

    m.chid_max = chi_d(m,m.dep_max,a)
    # println("chid: ", chid)
    # println("chi_0: ", m.chi_0)
    # println("chi_max: ", m.chid_max)
    # println("a_stable: ", m.a_stable)

    if chid > m.chid_max
        return m.dep_max
    elseif chid < -m.chid_max
        return -m.dep_max
    elseif chid > m.chi_0
        return (abs(chid - m.chi_0)/m.chi_1)^(1/m.chi_2)*m.a_stable
    elseif chid < -m.chi_0
        return -(abs(chid + m.chi_0)/m.chi_1)^(1/m.chi_2)*m.a_stable
    elseif !isnan(chid)
        println("!isnan case")
        return 0.0
    else
        return NaN
    end
end

function construct_AL!(m)
    m.AL = kron(sparse(I,m.na,m.na), sparse(I,m.nb, m.nb))
    m.AL = kron(m.yTtrans, m.AL)
    m.AL = kron(m.yPtrans, m.AL)
    m.AL = kron(m.yFtrans,m.AL)
    m.AL = kron(m.rhtrans, m.AL)
    m.AL -= I
    m.AL *= m.income_arrival
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
                    for rh in 1:m.nrh
                        # m.income_trans[a,yF,yP,yT,rh][a,:,:,:,:] = ones(m.nyF,m.nyP,m.nyT,m.nrh)./sum(ones(m.nyF,m.nyP,m.nyT,m.nrh))
                        
                        # Fixed
                        m.income_trans[a,yF,yP,yT,rh][a,:,:,:,:] = repeat(reshape(m.yFtrans[yF,:], (m.nyF,1,1,1)), 1, m.nyP, m.nyT, m.nrh)
                        # Persistent
                        m.income_trans[a,yF,yP,yT,rh][a,:,:,:,:] .*= repeat(reshape(m.yPtrans[yP,:], (1,m.nyP,1,1)), m.nyF, 1, m.nyT, m.nrh)
                        # Transitory
                        m.income_trans[a,yF,yP,yT,rh][a,:,:,:,:] .*= repeat(reshape(m.yTtrans[yT,:], (1,1,m.nyT,1)), m.nyF, m.nyP, 1, m.nrh)

                        # rho
                        m.income_trans[a,yF,yP,yT,rh][a,:,:,:,:] .*= repeat(reshape(m.rhtrans[rh,:], (1,1,1,m.nrh)), m.nyF, m.nyP, m.nyT, 1)
                    end
                end
            end
        end
    end
end

# Setup grids with power spacing
function power_grid(n,shape,lb,ub)
    if n == 1
        return [lb]
    end

    # Start with [0,1]
    grid = LinRange(0,1,n)
    # Powerspace
    grid = grid.^(1/shape)
    # Scale to bounds
    grid = lb .+ (ub- lb) .* grid

    return grid
    # # Apply to aagrid
    # m.aagrid = repeat(reshape(m.agrid, (m.na,1,1,1,1)), 1, m.nyF, m.nyP, m.nyT, m.nrh)

    # # Indicate that grids have been setup
    # m.grid_setup = true
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


# Returns normal approximation with optimal width
# Uses Roots.jl right now (may want to switch to my own solver)
function match_normal(n, mu, sigma, width_guess=0.5)
    width = find_zero(x -> discrete_normal(n,mu,sigma,x),width_guess)
    # temp, grid, dist = discrete_normal(n, mu, sigma, width)
    return discrete_normal(n, mu, sigma, width)
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



## 
#####################
# Test Kitchen
using Plots
##

m0 = ModelContinuous()
# m0.chi_2 = 0.15

a_test = 1.0
d_test = LinRange(-20,20,100)

plot(d_test, [chi(m0,d_t,a_test) for d_t in d_test])
# plot(d_test, [chi_d(m0,d_t,a_test) for d_t in d_test])
# plot(d_test, [chi_inv(m0,d_t,a_test) for d_t in d_test])
