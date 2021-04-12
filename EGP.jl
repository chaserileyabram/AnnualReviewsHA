include("ModelDiscrete.jl")

# Chase Abram
# For solving HA discrete-time models via endogenous grid point


function solve_EGP(m)

    # @unpack_Params p
    # unpack other structs also?

    # Guess cons
    con_old = m.r .* m.aagrid + m.yyFgrid + m.yyPgrid + m.yyTgrid
    # con_old = ones(size(m.aagrid))
    m.con = con_old

    # Initialize iteration trackers
    egp_iter = 0
    egp_diff = Inf

    println("Entering Loop...")
    # Run until convergence or out of iterations
    while egp_iter < m.egp_maxiter && egp_diff > m.egp_tol

        # Update each object in the param

        # u'(c)
        up(m)
        # println("it: ", egp_iter, ", var(m.up): ", var(m.up))

        # Expected marginal utility
        emuc(m)
        # println("it: ", egp_iter, ", var(m.emuc): ", var(m.emuc))

        # Marginal utility
        muc(m)
        # println("it: ", egp_iter, ", var(m.muc): ", var(m.muc))

        # Update forward looking consumption
        update_con_euler(m)
        # println("it: ", egp_iter, ", var(m.con_euler): ", var(m.con_euler))

        # Update forward looking assets
        update_a_euler(m)
        # println("it: ", egp_iter, ", var(m.a_euler): ", var(m.a_euler))

        # Interpolate/invert to get implied asset choice
        interp_a_tom(m)
        # println("it: ", egp_iter, ", var(m.a_tom): ", var(m.a_tom))

        # Implied consumption choice
        update_con(m)
        # println("it: ", egp_iter, ", var(m.con): ", var(m.con))

        
        # Check diff
        egp_diff = maximum(abs.(m.con - con_old))
        
        # Update
        con_old .= m.con
        egp_iter += 1

        # Add some print statements for progress here
        println("iter: ", egp_iter)
        println("    diff: ", egp_diff)
        
        # c_plot = plot(m.agrid, m.con[:,1,1,1,1],
        # xlabel = "a", ylabel = "c", title = "Consumption",
        # legend = false,)
        # display(c_plot)

        # a_plot = plot(m.agrid, m.a_tom[:,1,1,1,1],
        # xlabel = "a", ylabel = "a'", title = "Assets tomorrow",
        # legend = false,)
        # display(a_plot)


    end
end


# m0 = ModelDiscrete()
# setup_income(m0)
# # println("sum(m0.income_trans): ", sum(m0.income_trans))
# # println("size(m0.income_trans): ", size(m0.income_trans))
# # println("size(m0.income_trans[1,1,1,1,1]): ", size(m0.income_trans[1,1,1,1,1]))
# # println("con_old: ", m0.con)
# solve_EGP(m0)
# # println(m0.con)

# c_plot = plot(m0.agrid, m0.con[:,1,1,1,1],
# xlabel = "a", ylabel = "c", title = "Consumption",
# legend = false,)
# display(c_plot)

# a_plot = plot(m0.agrid, m0.a_tom[:,1,1,1,1], label = "a'",
# xlabel = "a", ylabel = "a'", title = "Assets tomorrow",
# legend = :bottomright)
# plot!(m0.agrid, m0.agrid, label = "45-deg")
# display(a_plot)

# println("a grid: ", m0.agrid)
# println("implied a'?: ", m0.a_tom[:,1,1,1,1])
