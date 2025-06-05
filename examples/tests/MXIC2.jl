using DoubleFloats
using LinearAlgebra
using ForwardDiff
using OrdinaryDiffEq

include("../pomin-types.jl")        # data type definitions
include("../pomin-io.jl")           # input/output routines
include("../idxer.jl")              # routines for index management
include("../integrators/rk4i.jl")   # Julia integrator functions
include("../HamPM.jl")              # Julia integrator functions
include("../integrators/intjul.jl")    # Julia integrator functions

#   Datatype for parameters
struct Parameters
    sym::Tuple{Bool,Real}     # Use the symplectic integrator?            The real parameter in the tuple is the ω parameter in the Tao map
    rkl::Tuple{Bool,Real}     # Use the rk4 integrator?                   The real parameter in the tuple is the Courant number
    jli::Tuple{Bool,Real}     # Use integrators in OrdinaryDiffEq.jl?     The real parameter in the tuple is the tolerance
    tspan::Tuple{Real,Real}   # Time span (initial_time,final_time)
    iter::Int                # Either the number of iterations, or maximum number of iterations, depending on integrator
end

m1  = Double64(0.01)
μ   = Double64(0)

p   = Double64(1/2)
b   = Double64(6e15)
d   = Double64(b*10^6)


"""
    dpfunc(p, b, d, m1, m2, G, c)

    Returns the analytical change in momentum for massless particle collision with:
    - p: magnitude of momentum
    - b: impact parameter
    - d: initial separation
    - m1: mass of particle 1
    - m2: mass of particle 2
    - G: gravitational constant (default 1)
    - c: speed of light (default 1)
"""
function dpfunc(p, b, d, m1, m2, G, c, tpfl=Float64)
    E1  = c*√((c*m1)^2 + p^2)
    E2  = c*√((c*m2)^2 + p^2)
    dp  = ( (tpfl(2)*G*(E1*E2)^2)/(b*p*(E1+E2)) )*
          ( tpfl(1) + ( tpfl(1)/E1^2 + tpfl(1)/E2^2 + tpfl(4)/(E1*E2) )*p^2 + p^4/(E1*E2)^2 )
    return dp
end

#-----------------------------------------------------------------------
function MXICgenML(p,b,dist;
                 G=Float64(1),c=Double64(1),
                 StartTime=Double64(0),MaxTimeStepFactor=Double64(1e5))

    n   = Double64(2)

    Dur = dist/c
    r1  = dist/Double64(2)
    r2  = dist/Double64(2)

    MaxTimeStep = Int(Dur/MaxTimeStepFactor)

    m1 = Double64(0)
    m2 = Double64(0)

    dp  = dpfunc(p, b, dist, m1, m2, G, c,Double64)

    qx1  = - r1
    qy1  = - b/Double64(2)
    qz1  = Double64(0.)

    qx2  = r2
    qy2  = b/Double64(2)
    qz2  = Double64(0.)

    px1  = p
    py1  = Double64(0.)
    pz1  = Double64(0.)

    px2  = -p
    py2  = Double64(0.)
    pz2  = Double64(0.)

    return ( dp , Parameters((false,Double64(1.)),
                             (true,Double64(0.005)),
                             (false,Double64(1e-6)),
                             (Double64(0.),StartTime+Dur),MaxTimeStep) ,
                    Double64[0.0;0.0] ,
                    [qx1;qy1;qz1;qx2;qy2;qz2;px1;py1;pz1;px2;py2;pz2] )
end #-------------------------------------------------------------------

# Test a range of nn values
for nn in [20, 40, 80, 100, 150]
    println("\nTesting with nn = ", nn)
    
    m1  = Double64(0.0);  # Massless particles
    μ   = Double64(π/4);
    db  = Double64(100);
    p   = Double64(1.0);  # Increased momentum for better numerical stability
    b   = Double64(10)*Double64("1.1108305558745590335689712446765042841434478759765625")^nn
    dst   = Double64(db*b);
    
    println("Impact parameter b = ", b)
    println("Initial separation dst = ", dst)
    
    mxicd=MXICgenML(p,b,dst;
                        G=Double64(1),c=Double64(1),
                        StartTime=Double64(0),MaxTimeStepFactor=Double64(1e4))  # Adjusted timestep factor
    
    # Wrap FHE in a form compatible with ODEProblem
    function wrapped_fhe!(du, u, p, t)
        du .= FHE(3, mxicd[3], u)
    end
    
    # Use moderate tolerances
    S = jlintegrator(mxicd[4], wrapped_fhe!, mxicd[2].tspan, nothing, Double64(1e-16), Double64(1e-16), AutoVern9(Rosenbrock23()))
    
    println("Analytical dp (mxicd[1]): ", mxicd[1])
    println("Numerical dp ((S[end]-S[1])[11]): ", (S[end]-S[1])[11])
    println("Percentage difference: ", Double64(100)*(abs((S[end]-S[1])[11])-mxicd[1])/mxicd[1], "%")
end

# hrkintegrator(3, 2, mxicd[4] , x -> pomin.dH_plus_MW(3, mxicd[3], x), δ, no_adapt, mxicd[2].tspan, mxicd[2].iter)
