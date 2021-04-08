using Parameters
using Random # only for testing

# Structure to hold all parameters and the setup functions of the problem

@with_kw struct Params

    # Discount factor
    beta0 = 0.98

    # Relative Risk Aversion
    gamma = 1.0

    # CRRA
    crra = 1 # 0 means Epstein-Zin

    na = 5

    
    nyF = 2
    nyP = 3
    nyT = 4

    nz = 3

    # income transitions (needs to be set up)
    income_trans = [ones(na,nyF,nyP,nyT,nz) for a in 1:na, yF in 1:nyF, yP in 1:nyP, yT in 1:nyT, z in 1:nz]

    # Expected MUC
    emuc = zeros(na,nyF,nyP,nyT,nz)


end

# Marginal utility
function up(p,c)
    if p.crra == 1
        return c.^(-p.gamma)
    else
        # Do Epstein-Zin
        return NaN
    end
end

# Expected MUC
function emuc(p,c)
    # Should add threading
    for a in 1:p.na
        for yF in 1:p.nyF
            for yP in 1:p.nyP
                for yT in 1:nyT
                    for z in 1:p.nz
                        p.emuc[a,yF,yP,yT,z] = sum(up(p,c) .* p.income_trans[a,yF,yP,yT,z])
                    end
                end
            end
        end
    end
end

p0 = Params()
# println(p0)

@unpack na,nyF,nyP,nyT,nz = p0
c0 = exp.(randn(na,nyF,nyP,nyT,nz))

up(p0, c0[1,1,1,1,1])
up(p0, c0)
emuc(p0,c0)
p0.emuc