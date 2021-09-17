#-----------------------------------------------------------------------
#   GENERAL INTEGRATOR FUNCTIONS
#-----------------------------------------------------------------------

#-----------------------------------------------------------------------

function IntODE( z0::RealVec , F::Function , tspan::Tuple{Real,Real} 
                , abstol::Real=1e-10 , reltol::Real=1e-10 
                , integrator=Tsit5() )

    f=(zx,p,t)->F(zx)

    prob = ODEProblem(f,z0,tspan,abstol,reltol)

    S = solve(prob,integrator)

    return (soln(S.t,S.u),S)
end

#-----------------------------------------------------------------------
