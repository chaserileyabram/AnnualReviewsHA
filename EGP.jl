


# Chase Abram
# For solving HA discrete-time models via endogenous grid point



function solve_EGP(p, grids, income)

    @unpack_Params
    # unpack other structs also?

    # Guess cons
    cons = r.*a_grid + y_grid

    egp_iter = 0
    egp_diff = Inf

    while egp_iter < egp_maxiter && egp_diff > egp_tol

        # Expected MUC
        

    end





