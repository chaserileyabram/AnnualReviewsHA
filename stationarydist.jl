# Chase Abram
# include("ModelDiscrete.jl")

# Just for now to test
include("EGP.jl")

# A = [0.9 0.1 0.0; 
# 0.0 0.6 0.4;
# 0.0 0.3 0.7]

# q = ones(3)
# q /= sum(q)

# it = 0

# while it < 1000
#     q = A'*q
#     q /= sum(q)
#     it += 1
# end
# println("q: ", q)

# Find stationary distribution over state space
function find_statdist(m)
    
    # Initialize
    statdist_iter = 0
    statdist_diff = Inf

    # Store old dist
    statdist_oldvec = m.statdistvec

    # Iterate until converged (or max iterations)
    while statdist_iter < m.statdist_maxiter && statdist_diff > m.statdist_tol
        
        # Update dist using adjoint of policies and transitions
        m.statdistvec = m.full_trans'*m.statdistvec

        # Normalize (corrects for flying off to infinities)
        m.statdistvec /= sum(m.statdistvec)

        # Get difference
        statdist_diff = maximum(abs.(m.statdistvec - statdist_oldvec))

        # Update
        statdist_oldvec = m.statdistvec
        statdist_iter += 1

        # Progress
        println("statdist_iter: ", statdist_iter)
        println("    diff: ", statdist_diff)
    end

    # Take vec back to grid
    m.statdist = reshape(m.statdistvec, (m.na, m.nyF, m.nyP, m.nyT, m.nbh))
end

# Find full transition matrix over states 
# Accounts for shocks and policy
function setup_full_trans(m)
    for a in 1:m.na
        for yF in 1:m.nyF
            for yP in 1:m.nyP
                for yT in 1:m.nyT
                    for bh in 1:m.nbh
                        # Fill in each row
                        m.full_trans[find_full_index(m,a,yF,yP,yT,bh), :] = vec(m.inter_trans[a,yF,yP,yT,bh])
                    end
                end
            end
        end
    end
    println("full_trans set up")
end

# Maps state index to column/row in full trans matrix
function find_full_index(m, a, yF, yP, yT, bh)
    @unpack_ModelDiscrete m
    return (bh-1)*na*nyF*nyP*nyT + (yT-1)*na*nyF*nyP + (yP-1)*na*nyF + (yF-1)*na + a
end

# Set up intermediate transmission matrix
# This is like the income transition, but accounts for policy also
function setup_inter_trans(m)
    # Zero out
    m.inter_trans *= 0.0

    # Add threading?
    for a in 1:m.na
        for yF in 1:m.nyF
            for yP in 1:m.nyP
                for yT in 1:m.nyT
                    for bh in 1:m.nbh

                        # Split implied a onto agrid
                        split_a_tom(m,m.a_tom[a,yF,yP,yT,bh])

                        # Use split to update transition
                        for a2 in 1:m.na
                            # Lower bound weight
                            if m.agrid[a2] == m0.lb
                                m.inter_trans[a,yF,yP,yT,bh][a2,:,:,:,:] += m.lbubmix .* m.income_trans[a,yF,yP,yT,bh][a,:,:,:,:]
                            end
                            # Upper bound weight
                            if m.agrid[a2] == m0.ub
                                m.inter_trans[a,yF,yP,yT,bh][a2,:,:,:,:] += (1 - m.lbubmix) .* m.income_trans[a,yF,yP,yT,bh][a,:,:,:,:]
                            end
                        end
                    end
                end
            end
        end
    end
    println("inter_trans set up")
end

# Returns the asset grid points that will recieve weight of policy ap
function split_a_tom(m, ap)

    # Different weights are used for on grid, 
    # above max, and above min, purely for diagnostics

    # Just use if on grid
    if ap in m.agrid
        m.lb = ap
        m.ub = ap
        m.lbubmix = 0.3
    # Use max if OOB
    elseif ap > maximum(m.agrid)
        m.lb = maximum(m.agrid)
        m.ub = maximum(m.agrid)
        m.lbubmix = 0.4
    # Use min if OOB
    elseif ap < minimum(m.agrid)
        m.lb = minimum(m.agrid)
        m.ub = minimum(m.agrid)
        m.lbubmix = 0.5
    # Else split weight
    else
        m.lb = maximum(m.agrid[ap .> m.agrid])
        m.ub = minimum(m.agrid[ap .< m.agrid])
        m.lbubmix = (m.ub - ap)/(m.ub - m.lb)
    end
end

function find_adist(m)
    m.adist = sum(m.statdist, dims=(2,3,4,5))[:,1,1,1,1]
end



m0 = ModelDiscrete(na = 50, nyF = 2, nyP = 2, 
nyT = 2, nbh = 2)
setup_power_grids(m0)
setup_income(m0)
solve_EGP(m0)

setup_inter_trans(m0)
setup_full_trans(m0)

find_statdist(m0)
find_adist(m0)

plot(m0.agrid, m0.adist)

# plot(m0.statdist)


# sum1 = sum(m0.income_trans, dims=2)
# println("sum 1: ", unique(sum(sum1)))
# sum2 = sum(m0.full_trans, dims=2)
# println("sum 2: ", unique(sum2))

# p1 = plot(sum1)
# display(p1)
# p2 = plot(sum2)
# display(p2)

# Need to start with income dist then 
# check each construction

# sum(m0.statdist, dims = (1,2))
# newstatdist = reshape(m0.statdist, (m0.na, m0.nyF, m0.nyP, m0.nyT, m0.nbh))

# newnewstatdist = sum(newstatdist, dims=(2,3,4,5))[:,1,1,1,1]
# p = plot(m0.agrid, newnewstatdist, title = "reasonable dist?")
# display(p)
# savefig(p, string(ENV["SLURM_ARRAY_TASK_ID"])* "sample_dist.png")


