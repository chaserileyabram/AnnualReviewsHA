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

# Get task id from slurm
if "SLURM_ARRAY_TASK_ID" in keys(ENV)
    task_id = parse(Int64, ENV["SLURM_ARRAY_TASK_ID"])
    println("SLURM_ARRAY_TASK_ID found with ", task_id)
else
    task_id = 9999
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

        setup_inter_trans(m)
        setup_full_trans(m)

        find_statdist(m)
        find_adist(m)

        return mean_wealth(m) - m.target_mean_wealth
    end

    println("lower: ", fitness(0.97))
    println("upper: ",fitness(0.99))
    # op = optimize(fitness, 0.92, 0.99)

    return fitness(0.98)

    # return find_zero(fitness, (0.9, 0.99), Bisection())
end

println("task_id (now creates calibrate fn): ", task_id)

m0 = ModelDiscrete(na = task_id)
println("m0 created with ", m0.na, " size agrid")
setup_power_grids(m0)
println("m0 grids setup")
setup_income(m0)
println("m0 income setup")

println("cal: ", calibrate(m0))

# println("m0 mean wealth target: ", m0.target_mean_wealth)
# println("calibrated mean wealth: ", mean_wealth(m0))


# f(x) = x^100
# op = optimize(f, -10, 10)

##




