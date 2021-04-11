include("Params.jl")

# Chase Abram
# For solving HA discrete-time models via endogenous grid point


function solve_EGP(p)

    # @unpack_Params p
    # unpack other structs also?

    # Guess cons
    con_old = p.r .* p.aagrid + p.yyFgrid + p.yyPgrid + p.yyTgrid
    # con_old = ones(size(p.aagrid))
    p.con = con_old

    # Initialize iteration trackers
    egp_iter = 0
    egp_diff = Inf

    println("Entering Loop...")
    # Run until convergence or out of iterations
    while egp_iter < p.egp_maxiter && egp_diff > p.egp_tol

        # Update each object in the param

        # u'(c)
        up(p)
        # println("it: ", egp_iter, ", var(p.up): ", var(p.up))

        # Expected marginal utility
        emuc(p)
        # println("it: ", egp_iter, ", var(p.emuc): ", var(p.emuc))

        # Marginal utility
        muc(p)
        # println("it: ", egp_iter, ", var(p.muc): ", var(p.muc))

        # Update forward looking consumption
        update_con_euler(p)
        # println("it: ", egp_iter, ", var(p.con_euler): ", var(p.con_euler))

        # Update forward looking assets
        update_a_euler(p)
        # println("it: ", egp_iter, ", var(p.a_euler): ", var(p.a_euler))

        # Interpolate/invert to get implied asset choice
        interp_a_tom(p)
        # println("it: ", egp_iter, ", var(p.a_tom): ", var(p.a_tom))

        # Implied consumption choice
        update_con(p)
        # println("it: ", egp_iter, ", var(p.con): ", var(p.con))

        
        # Check diff
        egp_diff = maximum(abs.(p.con - con_old))
        
        # Update
        con_old .= p.con
        egp_iter += 1

        # Add some print statements for progress here
        println("iter: ", egp_iter)
        println("    diff: ", egp_diff)
        
        # c_plot = plot(p.agrid, p.con[:,1,1,1,1],
        # xlabel = "a", ylabel = "c", title = "Consumption",
        # legend = false,)
        # display(c_plot)

        # a_plot = plot(p.agrid, p.a_tom[:,1,1,1,1],
        # xlabel = "a", ylabel = "a'", title = "Assets tomorrow",
        # legend = false,)
        # display(a_plot)


    end
end


p0 = Params()
setup_income(p0)
# println("sum(p0.income_trans): ", sum(p0.income_trans))
# println("size(p0.income_trans): ", size(p0.income_trans))
# println("size(p0.income_trans[1,1,1,1,1]): ", size(p0.income_trans[1,1,1,1,1]))
# println("con_old: ", p0.con)
solve_EGP(p0)
# println(p0.con)

c_plot = plot(p0.agrid, p0.con[:,1,1,1,1],
xlabel = "a", ylabel = "c", title = "Consumption",
legend = false,)
display(c_plot)

a_plot = plot(p0.agrid, p0.a_tom[:,1,1,1,1],
xlabel = "a", ylabel = "a'", title = "Assets tomorrow",
legend = false,)
display(a_plot)


