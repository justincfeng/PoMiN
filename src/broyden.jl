#-----------------------------------------------------------------------
#
#   BROYDEN ALGORITHM
#
#-----------------------------------------------------------------------

#-----------------------------------------------------------------------
"""
    JiSMU( ΔF::RealVec , Δx::RealVec , Ji::RealMtx )

The function `JiSMU` implements the Sherman-Morrison update formula, 
returning an updated value of the Jacobian matrix `Ji`, given the 
respective differences `ΔF` and `Δx` for the function ``F(x)`` and its 
argument ``x``.

"""
function JiSMU( ΔF::RealVec , Δx::RealVec , Ji::RealMtx )
    ΔxTJi = transpose(Δx)*Ji
    return Ji + ( ( Δx - Ji*ΔF )/( ΔxTJi*ΔF ) )*( ΔxTJi )
end     #---------------------------------------------------------------

#-----------------------------------------------------------------------
"""
    bsolve( F::Function , J::RealMtx , f0::RealVec , x0::RealVec ,
            nb::Int=24 )

The function `bsolve` implements the Broyden algorithm; in particular,
it finds the roots of the function `F(x)`, given an initial Jacobian 
matrix `J` and the initial guesses `f0` and `x0`. The parameter `nb` 
specifies the maximum number of iterations.

"""
function bsolve( F::Function , J::RealVec , f0::RealVec , x0::RealVec ,
                 nb::Int=24 )
    tpfl = typeof(x0[1])

    Ji = inv(J)
    println("Ji=",Ji)
    F0 = f0

    p = dot(F0,F0)

    ( b , k ) = ( true , 1 )

    if nb > 0
        Δx = - Ji*f0
        println("x0=",x0)
        println("Δx=",Δx)
        xi = x0 + Δx
        xs = xi
        if nb > 1
            println(stderr, "Broyden iteration # 1")
            Fi = F(xi)
            ΔF = Fi - F0
            F0 = Fi
            p  = norm(F0)
            println(stderr, "Norm of F = ", p)
            ps = p
            i  = 2
            while i <= nb && b
                Ji = JiSMU( ΔF , Δx , Ji )
                Δx = - Ji*F0
                xi = xi + Δx
                println(stderr, "Broyden iteration # ", i)
                Fi = F(xi)
                ( ΔF , F0 ) = ( Fi - F0 , Fi )
                pc = norm(F0)
                println(stderr, "Norm of F = ", pc)
                if pc<p && pc<ps
                    ( p , ps , xs ) = ( pc , pc , xi )
                elseif pc>p && k>=8
                    b = false
                else
                    p  = pc
                    k += 1
                end
                i += 1
            end
        end
    end

    return xs
end     #---------------------------------------------------------------
 