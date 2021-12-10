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

This function extracts particle positions from Z
"""
function Z2q( n::Int , d::Int , i::Int , Z::RealVec )
    tpfl = typeof(Z[1])
    q = zeros(tpfl,d)
    for j=1:d
        q[j] = Z[d*(i-1)+j]
    end
    return q
end

"""
    Z2p( n::Int , d::Int , i::Int , Z::RealVec )

This function extracts particle momenta from Z
"""
function Z2p( n::Int , d::Int , i::Int , Z::RealVec )
    tpfl = typeof(Z[1])
    p = zeros(tpfl,d)
    for j=1:d
        p[j] = Z[d*(i-1)+j+d*n]
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
