# These functions extract particle positions and momenta from Z

#function Z2q( n::Int , d::Int , i::Int , Z::RealVec )
#    tpfl = typeof(Z[1])
#    q = zeros(tpfl,d)
#    for j=1:d
#        q[j] = Z[d*(i-1)+j]
#    end
#    return q
#end

#function Z2p( n::Int , d::Int , i::Int , Z::RealVec )
#    tpfl = typeof(Z[1])
#    p = zeros(tpfl,d)
#    for j=1:d
#        p[j] = Z[d*(i-1+n)+j]
#    end
#    return p
#end

# These functions extract particle positions and momenta from Z
"""
    Z2q( n::Int , d::Int , i::Int , Z::RealVec )

This function extracts particle positions from Z for particle i
    n is the number of particles and d is the number of dimensions
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
end

"""
    Z2p( n::Int , d::Int , i::Int , Z::RealVec )

This function extracts particle momenta from Z for particle i
    n is the number of particles and d is the number of dimensions
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
end

function Part2m( Part::Particles )
    return Part.m
end

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
end
