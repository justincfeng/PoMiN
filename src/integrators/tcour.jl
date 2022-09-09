#-----------------------------------------------------------------------
#   SIMPLE ADAPTIVE TIMESTEPPING
#-----------------------------------------------------------------------

"""
    tcour( z::RealVec , C = 0.001 )

This function implements a simple adaptive timestepping function 
inspired by the Courant-Friedrichs-Condition. 
"""
function tcour( z::RealVec , C = 0.001 , d=3 )
    tpfl = typeof(z[1])
    n2 = length(z)

    if iseven(n2)
        n = Int(round(n2 / 2, digits=0))
        N = Int(round(n/d), digits=0)
        return 
    else
        return tpfl(0)
    end
end  # End tcour
