# Chase Abram
include("ModelDiscrete.jl")

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
    
    statdist_iter = 0
    statdist_diff = Inf

    statdist_old = m.statdist

    while statdist_iter < m.statdist_maxiter && statdist_diff > m.statdist_tol
        m.statdist = m.full_trans'*m.statdist

        m.statdist /= sum(m.statdist)

        statdist_diff = maximum(abs.(m.statdist - statdist_old))

        statdist_old = m.statdist
        statdist_iter += 1
        println("statdist_iter: ", statdist_iter)
        println("    diff: ", statdist_diff)
    end
end

# Find full transition matrix over states 
# Accounts for shocks and policy
function find_full_trans(m)
    for a in 1:m.na
        for yF in 1:m.nyF
            for yP in 1:m.nyP
                for yT in 1:m.nyT
                    for z in 1:m.nz
                        for a2 in 1:m.na
                            for yF2 in 1:m.nyF
                                for yP2 in 1:m.nyP
                                    for yT2 in 1:m.nyT
                                        for z2 in 1:m.nz
                                            m.full_trans[find_full_index(m,a, yF,yP,yT,z), find_full_index(m,a2, yF2,yP2,yT2,z2)] = m.inter_trans[a, yF,yP,yT,z][a2, yF2,yP2,yT2,z2]
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end

# Maps state index to column/row in full trans matrix
function find_full_index(m, a, yF, yP, yT, z)
    @unpack_ModelDiscrete m
    return (z-1)*na*nyF*nyP*nyT + (yT-1)*na*nyF*nyP + (yP-1)*na*nyF + (yF-1)*na + a
end

# Set up intermediate transmission matrix
# This is like the income transition, but accounts for policy also
function setup_inter_trans(m)
    m.inter_trans *= 0.0
    # Add threading?
    for a in 1:m.na
        for yF in 1:m.nyF
            for yP in 1:m.nyP
                for yT in 1:m.nyT
                    for z in 1:m.nz
                        # println("size(ones(yF,yP,yT,z)./sum(ones(yF,yP,yT,z))): ", size(ones(m.nyF,m.nyP,m.nyT,m.nz)./sum(ones(m.nyF,m.nyP,m.nyT,m.nz))))
                        # m.income_trans[a,yF,yP,yT,z][a2,:,:,:,:] = ones(m.nyF,m.nyP,m.nyT,m.nz)./sum(ones(m.nyF,m.nyP,m.nyT,m.nz))

                        # if m.a_tom[a, yF, yP, yT, z]

                        # end

                        for a2 in 1:m.na
                            m.acmix = (m.c - ap)/(m.c - m.a)

                            if a2 == split_a_tom(m,a)[1]
                                m.inter_trans[a,yF,yP,yT,z][a2,:,:,:,:] = m.acmix .* m.income_trans[a,yF,yP,yT,z][a,:,:,:,:]
                            elseif a2 == split_a_tom(m,a)[2]
                                m.inter_trans[a,yF,yP,yT,z][a2,:,:,:,:] = (1 - m.acmix) .* m.income_trans[a,yF,yP,yT,z][a,:,:,:,:]
                            end
                        end
                    end
                end
            end
        end
    end
end

# Returns the asset grid points that will recieve weight of policy ap
function split_a_tom(m, ap)
    if ap in m.agrid
        return [ap, ap]
    else
        if ap > maximum(agrid)
            m.lb = maximum(agrid)
            m.ub = maximum(agrid)
        elseif ap < minimum(agrid)
            m.lb = minimum(agrid)
            m.ub = minimum(agrid)
        else
            m.lb = maximum(m.agrid[ap .> m.agrid])
            m.ub = minimum(m.agrid[ap .< m.agrid])
        end
    end
        
    return [m.lb, m.ub]
end



m0 = ModelDiscrete()
# setup_income(m0)
# solve_EGP(m0)

# find_statdist(m0)
