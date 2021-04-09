using Params.jl

# Chase Abram
# For solving HA discrete-time models via endogenous grid point


function solve_EGP(p, grids, income)

    @unpack_Params p
    # unpack other structs also?

    # Guess cons
    con_old = r.*aa_grid + yyF_grid + yyP_grid + yyT_grid
    p.con = con_old

    # Initialize iteration trackers
    egp_iter = 0
    egp_diff = Inf

    # Run until convergence or out of iterations
    while egp_iter < egp_maxiter && egp_diff > egp_tol

        # Update each object in the param

        # u'(c)
        up(p)

        # [u'^{-1}]
        upinv(p)

        # Expected marginal utility
        emuc(p)

        # Marginal utility
        muc(p)

        # Update forward looking consumption
        update_con_euler(p)

        # Update forward looking assets
        update_a_euler(p)

        # Interpolate/invert to get implied asset choice
        interp_a_tom(p)

        # Implied consumption choice
        update_con(p)

        # Check diff
        egp_diff = maximum(abs.(p.con - con_old))

        # Update
        con_old = p.con
        egp_iter += 1

        # Add some print statements for progress here
        println("iter: ", egp_iter)
        println("    diff: ", egp_diff)
    end





