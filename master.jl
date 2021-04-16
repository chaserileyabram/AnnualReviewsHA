# Chase Abram

# Master
# Used for calibrating different model specs
# to data 
# using Pkg
# Pkg.add("Roots")
using Roots

using Optim
include("ModelDiscrete.jl")
include("EGP.jl")
include("stationarydist.jl")
include("statistics.jl")
include("mpc.jl")

##

function calibrate(m)
    
    setup_power_grids(m)
    println("m grids setup")
    setup_income(m)
    println("m income setup")
    
    calls = 0
    function fitness(b)
        call += 1
        println("fitness call: ", call)
        m.beta0 = b

        solve_EGP(m)

        setup_inter_trans(m)
        setup_full_trans(m)

        find_statdist(m)
        find_adist(m)

        return mean_wealth(m) - m.target_mean_wealth
    end

    println("lower: ", fitness(0.9))
    println("upper: ",fitness(0.99))
    # op = optimize(fitness, 0.92, 0.99)
    return find_zero(fitness, (0.9, 0.99), Bisection())
end


m0 = ModelDiscrete()
println("m0 created")

println("cal: ", calibrate(m0))

println("m0 mean wealth target: ", m0.target_mean_wealth)
println("calibrated mean wealth: ", mean_wealth(m0))


# f(x) = x^100
# op = optimize(f, -10, 10)




