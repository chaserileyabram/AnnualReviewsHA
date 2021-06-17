# HJB calculations

# Chase Abram

# Initial guess
function hjb_init!(m)
    m.V = uc(m, m.yyFgrid + m.yyPgrid + m.yyTgrid + m.rb .* m.bbgrid + m.ra .* m.aagrid)./(m.rho + m.zeta)
    # m.V = m.yyFgrid + m.yyPgrid + m.yyTgrid + m.rb .* m.bbgrid + m.ra .* m.aagrid
    m.Vvec = vec(m.V)
end

# Consumption
function consumption!(m)
    m.conF = ucp_inv(m,m.dVFb)
    m.conB = ucp_inv(m,m.dVBb)
    m.con0 = m.yyFgrid + m.yyPgrid + m.yyTgrid + m.rb.* m.bbgrid + m.ra .* m.aagrid
end

# Savings
function savings_c!(m)
    m.savcF = m.yyFgrid + m.yyPgrid + m.yyTgrid + m.rb.* m.bbgrid + m.ra .* m.aagrid - m.conF
    m.savcB = m.yyFgrid + m.yyPgrid + m.yyTgrid + m.rb.* m.bbgrid + m.ra .* m.aagrid - m.conB
end

# Hamiltonian
function hamiltonian_c!(m)
    m.HcF = uc(m,m.conF) + m.dVFb .* m.savcF
    m.HcB = uc(m,m.conB) + m.dVBb .* m.savcB
    m.Hc0 = uc(m, m.con0)
end

function upwind_con_indic!(m)
    # Need to account for NaNs on edges
    m.IcF = (m.conF .> 0) .& (m.savcF .> 0) .& ((m.savcB .>= 0) .| (((m.HcF .>= m.HcB) .| .!(m.conB .> 0)) .& ((m.HcF .>= m.Hc0) .| .!(m.con0 .> 0))))
    m.IcB = (m.conB .> 0) .& (m.savcB .< 0) .& ((m.savcF .<= 0) .| (((m.HcB .>= m.HcF) .| .!(m.conF .> 0)) .& ((m.HcB .>= m.Hc0) .| .!(m.con0 .> 0))))
    m.Ic0 = (m.con0 .> 0) .& .!(m.IcF .| m.IcB)
end

function upwind_con_apply!(m)
    m.con = m.IcF .* m.conF + m.IcB .* m.conB + m.Ic0 .* m.con0
    m.savc = m.IcF .* m.savcF + m.IcB .* m.savcB
end

# Deposits
function deposit!(m)
    # Loop over state space and fill FbFa, etc.
    for b in 1:m.nb
        for a in 1:m.na
            for yF in 1:m.nyF
                for yP in 1:m.nyP
                    for yT in 1:m.nyT
                        for rh in 1:m.nrh
                            # There should be NaNs in these
                            # Why are we getting depmaxes?

                            # if a == 1
                            #     println("m.dVFa[b,a,yF,yP,yT,rh]: ", m.dVFa[b,a,yF,yP,yT,rh])
                            #     println("m.dVBb[b,a,yF,yP,yT,rh]: ", m.dVBb[b,a,yF,yP,yT,rh])
                            # end

                            # println("chid: ", m.dVFa[b,a,yF,yP,yT,rh]/m.dVFb[b,a,yF,yP,yT,rh] - 1)
                            # println("chi_0: ", m.chi_0)
                            # println("a: ", m.aagrid[b,a,yF,yP,yT, rh])
                            m.depFbFa[b,a,yF,yP,yT, rh] = chi_inv(m,m.dVFa[b,a,yF,yP,yT,rh]/m.dVFb[b,a,yF,yP,yT,rh] - 1, m.aagrid[b,a,yF,yP,yT, rh])
                            m.depFbBa[b,a,yF,yP,yT, rh] = chi_inv(m,m.dVBa[b,a,yF,yP,yT,rh]/m.dVFb[b,a,yF,yP,yT,rh] - 1, m.aagrid[b,a,yF,yP,yT, rh])
                            m.depBbFa[b,a,yF,yP,yT, rh] = chi_inv(m,m.dVFa[b,a,yF,yP,yT,rh]/m.dVBb[b,a,yF,yP,yT,rh] - 1, m.aagrid[b,a,yF,yP,yT, rh])
                            m.depBbBa[b,a,yF,yP,yT, rh] = chi_inv(m,m.dVBa[b,a,yF,yP,yT,rh]/m.dVBb[b,a,yF,yP,yT,rh] - 1, m.aagrid[b,a,yF,yP,yT, rh])
                        end
                    end
                end
            end
        end
    end
end

function hamiltonian_d!(m)
    m.HdFF = m.dVFa .*m.depFbFa - m.dVFb .*(m.depFbFa + chi(m, m.depFbFa, m.aagrid))
    m.HdFB = m.dVBa .*m.depFbBa - m.dVFb .*(m.depFbBa + chi(m, m.depFbBa, m.aagrid))
    m.HdBF = m.dVFa .*m.depBbFa - m.dVBb .*(m.depBbFa + chi(m, m.depBbFa, m.aagrid))
    m.HdBB = m.dVBa .*m.depBbBa - m.dVBb .*(m.depBbBa + chi(m, m.depBbBa, m.aagrid))
end

function upwind_dep_indic!(m)

    # Initially set to check for validity
    # Need to correct so that NaNs don't mess up comparison

    # Never happens, but included for readability
    m.IdFF = (-m.depFbFa-chi(m,m.depFbFa, m.aagrid) .> 0) .& (m.depFbFa .> 0) #.& (m.HdFF .> 0)
    # These are likely
    m.IdFB = (-m.depFbBa-chi(m,m.depFbBa, m.aagrid) .> 0) .& (m.depFbBa .< 0) #.& (m.HdFB .> 0)
    m.IdBF = (-m.depBbFa-chi(m,m.depBbFa, m.aagrid) .< 0) .& (m.depBbFa .> 0) #.& (m.HdBF .> 0)
    # This one only happens if adjustment costs are high, but it is not chosen after convergence
    m.IdBB = (-m.depBbBa-chi(m,m.depBbBa, m.aagrid) .< 0) .& (m.depBbBa .< 0) #.& (m.HdBB .> 0)

    # Now actually set as indicators
    # (valid and others are either not valid, or worse Ham)
    m.IdFF = m.IdFF .& ((.!m.IdFB) .| (m.HdFF .>= m.HdFB)) .& ((.!m.IdFB) .| (m.HdFF .>= m.HdFB)) .& ((.!m.IdBB) .| (m.HdFF .>= m.HdBB))
    m.IdFB = m.IdFB .& ((.!m.IdFF) .| (m.HdFB .>= m.HdFF)) .& ((.!m.IdBF) .| (m.HdFB .>= m.HdBF)) .& ((.!m.IdBB) .| (m.HdFB .>= m.HdBB))
    m.IdBF = m.IdBF .& ((.!m.IdFF) .| (m.HdBF .>= m.HdFF)) .& ((.!m.IdFB) .| (m.HdBF .>= m.HdFB)) .& ((.!m.IdBB) .| (m.HdBF .>= m.HdBB))
    m.IdBB = m.IdBB .& ((.!m.IdFF) .| (m.HdBB .>= m.HdFF)) .& ((.!m.IdFB) .| (m.HdBB .>= m.HdFB)) .& ((.!m.IdBF) .| (m.HdBB .>= m.HdBF))
end

function upwind_dep_apply!(m)
    m.dep = m.IdFF .* m.depFbFa + m.IdFB .* m.depFbBa + m.IdBF .* m.depBbFa + m.IdBB .* m.depBbBa
    m.savd = -m.dep - chi(m, m.dep, m.aagrid)
end

function get_index(m,b,a,yF,yP,yT,rh)
    return (rh.-1).*m.nb.*m.na.*m.nyF.*m.nyP.*m.nyT .+ (yT.-1).*m.nb.*m.na.*m.nyF.*m.nyP .+ (yP.-1).*m.nb.*m.na.*m.nyF + (yF-1).*m.nb.*m.na .+ (a.-1).*m.nb .+ b
end

# Transitions used for solving HJB
function construct_A!(m)
    for b in 1:m.nb
        for a in 1:m.na
            for yF in 1:m.nyF
                for yP in 1:m.nyP
                    for yT in 1:m.nyT
                        for rh in 1:m.nrh
                            i = get_index(m,b,a,yF,yP,yT,rh)

                            m.A[i,i] = 0

                            # Liquid dissave
                            if b > 1
                                m.A[i, get_index(m,b-1,a,yF,yP,yT,rh)] = -(min(m.savc[b,a,yF,yP,yT,rh],0) + min(m.savd[b,a,yF,yP,yT,rh], 0))/(m.bbgrid[b,a,yF,yP,yT,rh] - m.bbgrid[b-1,a,yF,yP,yT,rh])
                                m.A[i,i] -= m.A[i, get_index(m,b-1,a,yF,yP,yT,rh)]
                            end

                            # Liquid save
                            if b < m.nb
                                m.A[i, get_index(m,b+1,a,yF,yP,yT,rh)] = (max(m.savc[b,a,yF,yP,yT,rh],0) + max(m.savd[b,a,yF,yP,yT,rh], 0))/(m.bbgrid[b+1,a,yF,yP,yT,rh] - m.bbgrid[b,a,yF,yP,yT,rh])
                                m.A[i,i] -= m.A[i, get_index(m,b+1,a,yF,yP,yT,rh)]
                            end

                            # Illiquid dissave
                            if a > 1
                                m.A[i,get_index(m,b,a-1,yF,yP,yT,rh)] = -(min(m.dep[b,a,yF,yP,yT,rh], 0) + min(m.ra*m.aagrid[b,a,yF,yP,yT,rh], 0))/(m.aagrid[b,a,yF,yP,yT,rh] - m.aagrid[b,a-1,yF,yP,yT,rh])
                                m.A[i,i] -= m.A[i,get_index(m,b,a-1,yF,yP,yT,rh)]
                            end

                            # Illiquid save
                            if a < m.na
                                m.A[i,get_index(m,b,a+1,yF,yP,yT,rh)] = (max(m.dep[b,a,yF,yP,yT,rh], 0) + max(m.ra*m.aagrid[b,a,yF,yP,yT,rh], 0))/(m.aagrid[b,a+1,yF,yP,yT,rh] - m.aagrid[b,a,yF,yP,yT,rh])
                                m.A[i,i] -= m.A[i,get_index(m,b,a+1,yF,yP,yT,rh)]
                            end
                        end
                    end
                end
            end
        end
    end
end

# Newton solver
function newton_solve(f,x; maxit = 100, tol = 1e-8)
    diff = Inf
    it = 0

    # x_new = NaN .* x
    while it < maxit && diff > tol
        x = x - ForwardDiff.jacobian(f,x) \ f(x)
        diff = maximum(abs.(f(x)))
        it += 1
        println("it: ", it, ", x: ", x)
    end
    return x
end

function update_Vvec!(m)
    # Need to set up AL and AD
    m.Vvec = ((m.rho + 1/m.hjb_delta).*sparse(I,m.full_trans_dim,m.full_trans_dim) - m.A - m.AL - m.AD) \ (uc(m,vec(m.con)) + 1/m.hjb_delta .* m.Vvec)
    m.V = reshape(m.Vvec, m.nb, m.na, m.nyF, m.nyP, m.nyT, m.nrh)
end

# Need to check corners

function iterate_hjb!(m)

    # Error
    # function error!(Vvec::AbstractArray)

        # m.Vvec = Vvec
        # m.V = reshape(m.Vvec, m.nb, m.na, m.nyF, m.nyP, m.nyT, m.nrh)

        # Update finite differences
        m.dVFb[1:end-1,:,:,:,:,:] = max.((m.V[2:end,:,:,:,:,:] - m.V[1:end-1,:,:,:,:,:])./(m.bbgrid[2:end,:,:,:,:,:] - m.bbgrid[1:end-1,:,:,:,:,:]), m.dVbmin)
        m.dVBb[2:end,:,:,:,:,:] = max.((m.V[2:end,:,:,:,:,:] - m.V[1:end-1,:,:,:,:,:])./(m.bbgrid[2:end,:,:,:,:,:] - m.bbgrid[1:end-1,:,:,:,:,:]), m.dVbmin)

        m.dVFa[:,1:end-1,:,:,:,:] = max.((m.V[:,2:end,:,:,:,:] - m.V[:,1:end-1,:,:,:,:])./(m.aagrid[:,2:end,:,:,:,:] - m.aagrid[:,1:end-1,:,:,:,:]), m.dVamin)
        m.dVBa[:,2:end,:,:,:,:] = max.((m.V[:,2:end,:,:,:,:] - m.V[:,1:end-1,:,:,:,:])./(m.aagrid[:,2:end,:,:,:,:] - m.aagrid[:,1:end-1,:,:,:,:]), m.dVamin)

        # Update consumption
        consumption!(m)
        # Update savings
        savings_c!(m)

        # Consumption Hamiltonian
        hamiltonian_c!(m)

        # Update upwind consumption indicators
        upwind_con_indic!(m)
        # Apply to cons-sav choice
        upwind_con_apply!(m)

        # Update deposits (all 4)
        deposit!(m)

        # Deposit Hamiltonians
        hamiltonian_d!(m)

        # Update upwind deposit indicators
        upwind_dep_indic!(m)
        # Apply to deposit choice
        upwind_dep_apply!(m)

        # Build implied A matrix
        construct_A!(m)
        # println("m.A: ", m.A)

        update_Vvec!(m)

        # Return error in HJB
        # return m.rho*m.Vvec - vec(uc(m,m.con)) - (m.A + m.AL + m.AD)*m.Vvec

        # f!(v) = m.rho*v - vec(uc(m,m.con)) - (m.A + m.AL + m.AD)*m.Vvec

        # m.Vvec = nlsolve()
        # m.V = reshape(m.V,m.nb, m.na, m.nyF, m.nyP, m.nyT, m.nrh)
    # end

    # return error!(m.Vvec)
    # m.Vec = newton_solve(error!, m.Vvec)
    # m.Vec = nlsolve(error!, m.Vvec, show_trace = true)
end

function solve_hjb!(m)
    runavg_diff = 1
    diff = Inf
    it = 0
    old_V = m.V
    while it < m.hjb_maxiter && runavg_diff > m.hjb_tol
        iterate_hjb!(m)
        diff = maximum(abs.(m.V - old_V))
        runavg_diff = .1 *runavg_diff + 0.9 * diff
        println("it: ", it, ", diff(raw): ", diff)
        println("it: ", it, ", runavg_diff(raw): ", runavg_diff)

        # diff = max(diff, maximum(abs.(m.rho*m.Vvec - vec(uc(m,m.con)) - (m.A + m.AL + m.AD)*m.Vvec)))
        # println("it: ", it, ", diff(PDE): ", diff)

        old_V = m.V
        it += 1
    end

    println("error: ", maximum(abs.(m.rho*m.Vvec - vec(uc(m,m.con)) - (m.A + m.AL)*m.Vvec)))
end

##
#################
# Test Kitchen

m0 = ModelContinuous()
println("m0 initialized")

##
hjb_init!(m0)
construct_AL!(m0)

# # println("before: ", m0.depFbFa)
# hjb_update!(m0)
# # println("after: ", m0.depFbFa)

# minimum(m0.IdFF + m0.IdFB + m0.IdBF + m0.IdBB)
# sum(m0.IdFF + m0.IdFB + m0.IdBF + m0.IdBB)

# sum(m0.dep)/m0.full_trans_dim

# f(x) = [(x[1] - 1)^6, (x[2] - 2)^2, (x[3] + 3)^2, (x[4]+ 300)^2]

# # nlsolve(f,[0.0, 0.0, 0.0])
# newton_solve(f, [0,0,0,0])

##

# iterate_hjb!(m0)
solve_hjb!(m0)

function count_nan(a)
    return sum(isnan.(a))
end

##
# surface(m0.agrid, m0.bgrid, m0.con[:,:,1,1,1,1])
plot(m0.bgrid[1:100], m0.con[1:100,1,1,1,1,1])
# plot(m0.agrid[1:11], m0.con[1,1:11,1,1,1,1])




