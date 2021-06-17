# Chase Abram

# Master
# Used for calibrating different model specs
# to data 
# using Pkg
# Pkg.add("Roots")
using Roots

using Optim

include("roots_HA.jl")

include("ModelDiscrete.jl")
include("EGP.jl")
include("stationarydist.jl")
include("statistics.jl")
include("mpc.jl")

##

# Get task id from slurm
if "SLURM_ARRAY_TASK_ID" in keys(ENV)
    task_id = parse(Int64, ENV["SLURM_ARRAY_TASK_ID"])
    println("SLURM_ARRAY_TASK_ID found with ", task_id)
else
    task_id = 50
end


function calibrate(m)
    
    setup_power_grids(m)
    println("m grids setup")
    setup_income(m)
    println("m income setup")
    
    call = 0
    function fitness(b)
        call += 1
        println("fitness call: ", call)
        m.beta0 = b

        solve_EGP(m)
        # println("EGP solved")
        # println("m.up: ", m.up)
        # println("m.emuc: ", m.emuc)
        # println("m.muc: ", m.muc)
        # println("m.con_euler: ", m.con_euler)
        # println("m.a_euler: ", m.a_euler)

        # # u'(c)
        # up(m)

        # # Expected marginal utility
        # emuc(m)

        # # Marginal utility
        # muc(m)

        # # Update forward looking consumption
        # update_con_euler(m)

        # # Update forward looking assets
        # update_a_euler(m)

        # # Interpolate/invert to get implied asset choice
        # interp_a_tom(m)

        # # Implied consumption choice
        # update_con(m)

        # println("m.a_tom: ", m.a_tom)

        setup_inter_trans(m)
        # println("inter_trans setup")
        setup_full_trans(m)
        # println("full_trans setup")

        find_statdist(m)
        # println("statdist found")
        # find_adist(m)
        # println("adist found")

        # return 99999

        return mean_wealth(m) - m.target_mean_wealth
    end

    lower = 0.95
    upper = 0.99
    println("lower: ", fitness(lower))
    println("upper: ",fitness(upper))
    # op = optimize(fitness, 0.92, 0.99)

    # return fitness(0.98)

    # bs = collect(LinRange(0.96,0.999, 100))
    # p = plot(bs, fitness.(bs),
    # xlabel = "beta", ylabel = "mean_w - target_w",
    # title = "Dev from target")
    # display(p)

    # return NaN

    # Should I build my own zero finder to optimize to the problem at hand?
    # return find_zero(fitness, (0.97, 0.99), 
    # Roots.Brent(); atol = 1e-8, rtol = 1e-8)

    final_beta = bisect(fitness, lower, upper, 0.5, 0.5;
    maxit = 15, ftol = 1e-6)

    fitness(final_beta)

    return final_beta

    # return newton(fitness, lower, upper; 
    # maxit = 20, ftol = 1e-6,)
end

println("task_id (now creates calibrate fn): ", task_id)

m0 = ModelDiscrete(na = task_id, a_max = 200)
println("m0 created with ", m0.na, " size agrid")
setup_power_grids(m0)
println("m0 grids setup")
setup_income(m0)
println("m0 income setup")
println("m0.agrid: ", m0.agrid)
println("m0.yFgrid: ", m0.yFgrid)
println("m0.yPgrid: ", m0.yPgrid)
println("m0.yTgrid: ", m0.yTgrid)


cal0 = calibrate(m0)
# m0.beta0 = cal0

println("cal: ", cal0)

println("m0 mean wealth target: ", m0.target_mean_wealth)
println("calibrated mean wealth: ", mean_wealth(m0))


# f(x) = x^100
# op = optimize(f, -10, 10)

##

# Need to plot mean_wealth function to see how best to find root - Done?

# Need to plot MPC as fn of beta. How sensitive?





