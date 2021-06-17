include("ModelDiscrete.jl")

# Chase Abram
# For solving HA discrete-time models via endogenous grid point


function solve_EGP(m)

    # Guess cons
    con_old = m.r .* m.aagrid + m.yyFgrid + m.yyPgrid + m.yyTgrid
    m.con = con_old

    # Initialize iteration trackers
    egp_iter = 0
    egp_diff = Inf

    # Should probably add "silent" option to model
    # println("Entering Loop...")

    # Run until convergence or out of iterations
    while egp_iter < m.egp_maxiter && egp_diff > m.egp_tol

        # Update each object in the param

        # u'(c)
        up(m)

        # Expected marginal utility
        emuc(m)

        # Marginal utility
        muc(m)

        # Update forward looking consumption
        update_con_euler(m)

        # Update forward looking assets
        update_a_euler(m)

        # Interpolate/invert to get implied asset choice
        interp_a_tom(m)

        # Implied consumption choice
        update_con(m)
        
        # Check diff
        egp_diff = maximum(abs.(m.con - con_old))
        
        # Update
        con_old .= m.con
        egp_iter += 1

        # Add some print statements for progress here
        # println("egp_iter: ", egp_iter)
        # println("    egp_diff: ", egp_diff)
        
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
# setup_power_grids(m0)
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
