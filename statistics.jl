# Chase Abram
# For computing statistics

include("ModelDiscrete.jl")
include("stationarydist.jl")
using Statistics
# using Pkg
# Pkg.add("QuadGK")
using QuadGK

# Mean over all types of income
function mean_income(m)
    println("mean_income: ", sum(m.statdist .* (m.yyFgrid .* m.yyPgrid .* m.yyTgrid)))
    return sum(m.statdist .* (m.yyFgrid .* m.yyPgrid .* m.yyTgrid))
end

# Variance in log of income (includes all types)
function var_log_gross_income(m)
    println("var_log_gross_income: ", sum(m.statdist .* log.(m.yyFgrid .* m.yyPgrid .* m.yyTgrid).^2) - sum(m.statdist .* log.(m.yyFgrid .* m.yyPgrid .* m.yyTgrid))^2)
    return sum(m.statdist .* log.(m.yyFgrid .* m.yyPgrid .* m.yyTgrid).^2) - sum(m.statdist .* log.(m.yyFgrid .* m.yyPgrid .* m.yyTgrid))^2
end

# Mean wealth
function mean_wealth(m)
    # println("mean_wealth: ", sum(m.statdist .* m.aagrid))
    return sum(m.statdist .* m.aagrid)
end

# Prob. wealth <= a
function wealth_lesseq(m,a)
    # println("wealth less than or equal to ", a, ": ", sum(m.statdist[m.aagrid .<= a]))
    return sum(m.statdist[m.aagrid .<= a])
end

# Prob. sav <= a
function a_tom_lesseq(m,a)
    println("savings less than or equal to ", a, ": ", sum(m.statdist[m.a_tom .<= a]))
    return sum(m.statdist[m.a_tom .<= a])
end

# Prob. wealth <= f*income
function wealth_lesseq_fi(m,f)
    println("wealth less than or equal to ", f, "*income: ", sum(m.statdist[m.aagrid .<= (f .* (m.yyFgrid .* m.yyPgrid .* m.yyTgrid))]))
    return sum(m.statdist[m.aagrid .<= (f .* (m.yyFgrid .* m.yyPgrid .* m.yyTgrid))])
end

# Wealth level at quantile q
function wealth_quantile(m,q)
    if q <= wealth_lesseq(m,m.agrid[1])
        wbc = wealth_lesseq(m,m.agrid[1])
        # println("Requested ", q,", but ", wbc, " at borrowing constraint")
        return 0
    end

    m.tmp_itp = interpolate((cumsum(m.adist),) , m.agrid, Gridded(Linear()))
    # println("wealth at quantile ", q, ": ", m.tmp_itp(q))
    return extrapolate(m.tmp_itp, Line())(q)
end

# Wealth held by bottom q
function wealth_held_bottom(m, q)
    # println("Bottom ", q, " holds wealth: ", sum((m.statdist .* m.aagrid)[m.aagrid .<= wealth_quantile(m,q)])/sum(m.statdist .* m.aagrid))
    return sum((m.statdist .* m.aagrid)[m.aagrid .<= wealth_quantile(m,q)])/sum(m.statdist .* m.aagrid)
end

# Wealth held by top q
function wealth_held_top(m, q)
    println("Top ", q, " holds wealth: ", 1 - wealth_held_bottom(m,1-q))
    return 1 - wealth_held_bottom(m,1-q)
end

# Gini coefficient
# 0 -> Equality
# 1 -> Inequality
function gini(m)
    # println("gini: ", 1 - sum(cumsum(sum(m.statdist, dims=(2,3,4,5))[:,1,1,1,1])))
    # return 1 - 2*sum(cumsum(sum(m.statdist, dims=(2,3,4,5))[:,1,1,1,1]))

    lorenz(z) = wealth_held_bottom(m,z)

    ps = LinRange(0,1,100)

    pl = plot(ps, lorenz.(ps), title = "Lorenz", 
    label = "The Sad Reality", legend = :topleft)
    plot!(ps, ps, label = "Utopian Equality")
    display(pl)
    println("gini: ", 1 - 2*quadgk(lorenz, 0, 1, rtol=1e-6)[1])
    return 1 - 2*quadgk(lorenz, 0, 1, rtol=1e-6)[1]
end


m0 = ModelDiscrete(beta0 = 0.5)
setup_power_grids(m0)
setup_income(m0)

solve_EGP(m0)

setup_inter_trans(m0)
setup_full_trans(m0)

find_statdist(m0)
find_adist(m0)


mean_income(m0)

var_log_gross_income(m0)

mean_wealth(m0)

wealth_lesseq(m0,0.0)
wealth_lesseq(m0,1.0)

a_tom_lesseq(m0, 0.0)
a_tom_lesseq(m0, 1.0)

wealth_lesseq_fi(m0,0)
wealth_lesseq_fi(m0,1/6)
wealth_lesseq_fi(m0,1/12)

wealth_quantile(m0, 0.0)
wealth_quantile(m0, 0.1)
wealth_quantile(m0, 0.25)
wealth_quantile(m0, 0.5)
wealth_quantile(m0, 0.9)
wealth_quantile(m0, 0.99)
wealth_quantile(m0, 0.999)

wealth_held_top(m0, 1)
wealth_held_top(m0, 0.5)
wealth_held_top(m0, 0.1)
wealth_held_top(m0, 0.01)

gini(m0)

