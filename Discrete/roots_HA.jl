# Chase Abram

# Home-grown root finders

function bisect(f, xl, xh, lw, uw; 
    maxit = 100, ftol = 1e-8)

    # Need to keep best

    fl = f(xl)
    fh = f(xh)

    xbest = 0
    fbest = Inf
    
    if abs(fl) < abs(fh)
        xbest = xl
        fbest = fl
    else
        xbest = xh
        fbest = fh
    end

    sl = sign(fl)
    sh = sign(fh)
    
    # Check if bracketing
    if sl == sh
        println("Improper initial bracket in bisect")
        return NaN
    end

    it = 0
    fdiff = Inf

    w = 1/2
    xm = 0
    fm = 0
    sm = 0

    while it < maxit && fdiff > ftol

        # Adjust weighting
        # if it % 10 == 1
        #     w = lw
        # elseif it % 10 == 2
        #     w = uw
        # else
        #     w = min(max(lw, 1/(exp(-fm) + 1)), uw)
        # end
        # w = 1/(exp(-fm) + 1)

        w = min(max(lw, 1/(exp(-fm) + 1)), uw)

        xm = w*xl + (1-w)*xh #+ (it % 10 == 1)*1e-6 - (it % 10 ==2)*1e-6
        fm = f(xm)
        sm = sign(fm)
        fdiff = abs(fm)

        if fdiff < abs(fbest)
            xbest = xm
            fbest = fm
            println("xbest updated")
            println("    xbest: ", xbest)
            println("    fbest: ", fbest)
        end

        println("it: ", it)
        println("    w: ", w)
        # println("    fm: ", fm)
        # println("    fdiff: ", fdiff)
        if sm == sl
            xl = xm
            fl = fm
            sl = sm
        else
            xh = xm
            fh = fm
            sh = sm
        end
        println("    fl: ", fl)
        println("    fm: ", fm)
        println("    fh: ", fh)


        it += 1
    end

    if it == maxit
        println("Failed, but within ", fbest)
    end

    println("final xbest: ", xbest)

    return xbest

end


# This one sucks
function secant(f,xl, xh; maxit = 100, ftol = 1e-8)
    
    fl = f(xl)
    fh = f(xh)

    sl = sign(fl)
    sh = sign(fh)
    
    # Check if bracketing
    if sl == sh
        println("Improper initial bracket in bisect")
        return NaN
    end

    it = 0
    fdiff = Inf

    xm = 0
    fm = 0
    sm = 0

    while it < maxit && fdiff > ftol

        xm = xl - fl*(xh-xl)/(fh-fl) #+ (it % 10 == 1)*1e-6 - (it % 10 ==2)*1e-6
        fm = f(xm)
        sm = sign(fm)
        fdiff = abs(fm)

        println("it: ", it)
        # println("    w: ", w)
        # println("    fm: ", fm)
        # println("    fdiff: ", fdiff)
        if sm == sl
            xl = xm
            fl = fm
            sl = sm
        else
            xh = xm
            fh = fm
            sh = sm
        end
        println("    fl: ", fl)
        println("    fm: ", fm)
        println("    fh: ", fh)


        it += 1
    end

    return xm
end


function newton(f, xl, xh; 
    maxit = 100, ftol = 1e-8)

    fl = f(xl)
    fh = f(xh)

    fbest = Inf

    if abs(fl) < abs(fh)
        xbest = xl
        fbest = abs(fl)
    else
        xbest = xh
        fbest = abs(fh)
    end

    sl = sign(fl)
    sh = sign(fh)
    
    # Check if bracketing
    if sl == sh
        println("Improper initial bracket in newton")
        return NaN
    end

    it = 0
    fdiff = Inf

    xm = 0
    fm = 0
    sm = 0
    Dfm = 0

    w = 1/2

    while it < maxit && fdiff > ftol

        # First bisect
        xm = (xl + xh)/2

        fm = f(xm)
        sm = sign(fm)

        if abs(fm) < fbest
            xbest = xm
            fbest = abs(fm)
        end

        if sm == sl
            Dfm = (fm - fl)/(xm - xl)
            xl = xm
            fl = fm
            sl = sm
        else
            Dfm = (fh - fm)/(xh - xm)
            xh = xm
            fh = fm
            sh = sm
        end

        # Try to Newton out
        xnew = xm - fm/Dfm

        if (xnew > xl && xnew < xh)
            println("newton success")
            xm = xnew
            fm = f(xm)
            sm = sign(fm)
        elseif xnew <= xl
            println("shot below xl")
            w = 1 - abs(2/(exp(-Dfm)+1)-1)
            xm = w*xl + (1-w)*xh
            fm = f(xm)
            sm = sign(fm)
        else
            println("shot above xh")
            w = abs(2/(exp(-Dfm)+1)-1)
            xm = w*xl + (1-w)*xh
            fm = f(xm)
            sm = sign(fm)
        end

        # else use Dfm to inform where to choose xm?

        # if it > maxit/2
        #     println("newton risky")
        #     xm = xnew
        #     fm = f(xm)
        #     sm = sign(fm)
        # end

        if abs(fm) < abs(fbest)
            xbest = xm
            fbest = abs(fm)
        end

        if sm == sl
            xl = xm
            fl = fm
            sl = sm
        else
            xh = xm
            fh = fm
            sh = sm
        end

        fdiff = abs(fm)

        println("it: ", it)
        println("    fl: ", fl)
        println("    fm: ", fm)
        println("    fh: ", fh)
        println("    Dfm: ", Dfm)
        it += 1
    end

    if it == maxit
        println("Failed, but within ", fbest)
    end

    return xbest

end



f(x) = (exp(x) - exp(5))

println("bisect: ", bisect(f, 0, 11, 0.25, 0.75; ftol = 1e-9))
# f(bisect(f, 0, 11; ftol = 1e-9))
# println("secant: ", secant(f,0, 11; ftol = 1e-8))
# println("newton: ", newton(f, 0, 11; ftol = 1e-8))
