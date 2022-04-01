#-----------------------------------------------------------------------
#   GENERAL INTEGRATOR FUNCTIONS
#-----------------------------------------------------------------------

#-----------------------------------------------------------------------

"""
    hointegratorfull( z0::RealVec , F::Function 
                        , tspan::Tuple{Real,Real} , abstol::Real=1e-10 
                        , reltol::Real=1e-10 , integrator=Tsit5() )

This function performs the integration using the integrators provided
in the Julia library 'OrdinaryDiffEq.jl', The default integrator is 
'Tsit5()'.
"""
function hointegratorfull( z0::RealVec , F::Function 
                        , tspan::Tuple{Real,Real} , abstol::Real=1e-10 
                        , reltol::Real=1e-10 , integrator=Tsit5() )

    f=(zx,p,t)->F(zx)

    prob = ODEProblem(f,z0,tspan,abstol,reltol)

    S = solve(prob,integrator)

    return (soln(S.t,S.u),S)
end

function hointegrator( z0::RealVec , F::Function 
                    , tspan::Tuple{Real,Real} , abstol::Real=1e-10 
                    , reltol::Real=1e-10 , integrator=Tsit5() )

    f=(zx,p,t)->F(zx)

    prob = ODEProblem(f,z0,tspan,abstol,reltol)

    S = solve(prob,integrator)

    return soln(S.t,S.u)
end

#-----------------------------------------------------------------------
