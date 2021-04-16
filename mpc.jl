# Chase Abram

# MPCS

# MPC of injection of delta (into cash-on-hand)
function mpc(m, delta)

    if delta == 0
        println("Cannot compute infinitesimal MPC")
        return NaN
    end

    mpc_out = 0
    
    for yF in 1:m.nyF
        for yP in 1:m.nyP
            for yT in 1:m.nyT
                for bh in 1:m.nbh
                    # Interpolate: x -> c
                    m.tmp_itp = interpolate((m.R .* m.aagrid[:,yF,yP,yT,bh] + m.yyFgrid[:,yF,yP,yT,bh] .* m.yyPgrid[:,yF,yP,yT,bh] .* m.yyTgrid[:,yF,yP,yT,bh],), m.con[:,yF,yP,yT,bh], Gridded(Linear()))
                    # Extrapolate for those at top and bottom
                    m.tmp_itp = extrapolate(m.tmp_itp, Line())
                    # Add to total (with dist weighting)
                    mpc_out += sum(m.statdist[:,yF,yP,yT,bh] .* (m.tmp_itp.(m.agrid .+ delta) - m.tmp_itp.(m.agrid))./delta)
                end
            end
        end
    end
    println("mpc_out of ", delta, ": ", mpc_out)
    return mpc_out

    # m.tmp_itp = interpolate((m.R .* m.aagrid + m.yyFgrid .* m.yyPgrid .* m.yyTgrid,), (m.con,), Gridded(Linear()))
    # m.tmp_itp = extrapolate(m.tmp_itp, Line())
    # return sum(m.statdist .*(m.tmp_itp.(m.R .* m.aaagrid .+ delta, m.yyFgrid, m.yyPgrid, m.yyTgrid, m.bbhgrid) - m.tmp_itp.(m.R .* m.aaagrid, m.yyFgrid, m.yyPgrid, m.yyTgrid, m.bbhgrid)))
end

# mpc(m0,0.1)
# mpc(m0,-0.1)

# mpc(m0,1.0)
# mpc(m0,-1.0)

# mpc(m0,10.0)
# mpc(m0,-10.0)

# deltas = collect(LinRange(-1,1,50))
# filter!(z -> z != 0, deltas)
# plot(deltas, [mpc(m0, d) for d in deltas], 
# title = "E[MPC(delta)]", xlabel = "delta", 
# ylabel = "E[m]", legend = false)