#-----------------------------------------------------------------------
# These functions extract particle positions and momenta from Z
#-----------------------------------------------------------------------

"""
    Z2q( n::Int , d::Int , i::Int , Z::RealVec )

This function extracts particle positions from Z. The parameter `n` is
the number of particles, `d` represents the spatial dimension, and `i`
is the particle label.
"""
function Z2q( n::Int , d::Int , i::Int , Z::RealVec )
    if i > n
        error("In Z2q: particle number i=" * string(i) * " exceeds number of particles n=" * string(n))
    end
    tpfl = typeof(Z[1])
    q = zeros(tpfl,d)
    for j=1:d
        index = d * (i - 1) + j
        q[j] = Z[index]
    end
    return q
end #-------------------------------------------------------------------

"""
    Z2p( n::Int , d::Int , i::Int , Z::RealVec )

This function extracts particle momenta from Z. The parameter `n` is
the number of particles, `d` represents the spatial dimension, and `i`
is the particle label.
"""
function Z2p( n::Int , d::Int , i::Int , Z::RealVec )
    if i > n
        error("In Z2p: particle number i=" * string(i) * " exceeds number of particles n=" * string(n))
    end
    tpfl = typeof(Z[1])
    p = zeros(tpfl,d)
    for j=1:d
        index = d * (i - 1) + j + d * n
        p[j] = Z[index]
    end
    return p
end #-------------------------------------------------------------------

"""
   Part2m( Part::Particles )

This function extracts particle masses from the `Particles` datatype.
"""
function Part2m( Part::Particles )
    return Part.m
end #-------------------------------------------------------------------

"""
   Part2Z( Part::Particles )

This function extracts the phase space vector from the `Particles`
datatype.
"""
function Part2Z( Part::Particles )
    tpfl = typeof(Part.q[1][1])
    n = length(Part.m)
    d = length(Part.q[1])

    if n==length(Part.p) && d==length(Part.p[1])
        Z = zeros(tpfl,2*n*d)
            for i=1:n
                for j=1:d
                    Z[d*(i-1)+j]   = Part.q[i][j]
                    Z[d*(i-1+n)+j] = Part.p[i][j]
                end
            end
        return Z
    else
        print("Inputs have inconsistent dimensionality \n")
    end
end #-------------------------------------------------------------------

"""
    CombineParticles2Z( Part1::Particles, Part2::Particles )

This function combines two `Particles` datatypes into a single phase space vector.
"""
function CombineParticles2Z( Part1::Particles, Part2::Particles )
    tpfl = typeof(Part1.q[1][1])
    n1 = length(Part1.m)
    n2 = length(Part2.m)
    d = length(Part1.q[1])

    if n1==length(Part1.p) && d==length(Part1.p[1]) && n2==length(Part2.p) && d==length(Part2.p[1])
        Z = zeros(tpfl,2*(n1+n2)*d)
            for i=1:n1
                for j=1:d
                    Z[d*(i-1)+j]   = Part1.q[i][j]
                    Z[d*(i-1+n1)+j] = Part1.p[i][j]
                end
            end
            for i=1:n2
                for j=1:d
                    Z[d*(i-1+n1)+j]   = Part2.q[i][j]
                    Z[d*(i-1+n1+n2)+j] = Part2.p[i][j]
                end
            end
        return Z
    else
        print("Inputs have inconsistent dimensionality \n")
    end
end #-------------------------------------------------------------------

#-----------------------------------------------------------------------
#   SYMPLECTIC OPERATOR
#-----------------------------------------------------------------------
"""
    Jsympl( Zarg::RealVec )

Symplectic operator. Maps output of dH to time derivative of phase space
variables.
"""
function Jsympl( Zarg::RealVec )
    tpfl=typeof(Zarg[1])
    n2 = length(Zarg)
    Z = zeros(tpfl,n2)

    if iseven(n2)
        n = Int(round(n2 / 2, digits=0))
        for i=1:n
            Z[i]    = Zarg[n+i] 
            Z[n+i]  = - Zarg[i] 
        end
        return Z
    else
        return Z
    end
end #-------------------------------------------------------------------