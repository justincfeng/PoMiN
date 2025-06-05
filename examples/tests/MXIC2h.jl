using ArbNumerics
using LinearAlgebra
using ForwardDiff
using OrdinaryDiffEq

setworkingprecision(ArbFloat, 128) # Set high precision

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

m1  = ArbFloat("0.01")
μ   = ArbFloat(0)

p   = ArbFloat("0.5")
b   = ArbFloat("6e15")
d   = b*ArbFloat("1e6")


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
                 G=ArbFloat(1),c=ArbFloat(1),
                 StartTime=ArbFloat(0),MaxTimeStepFactor=ArbFloat("1e5"))

    n   = ArbFloat(2)

    Dur = dist/c
    r1  = dist/ArbFloat(2)
    r2  = dist/ArbFloat(2)

    # Cap MaxTimeStep at Int64 max to avoid overflow
    MaxTimeStep = min(Int(1e6), Int(floor(min(Dur/MaxTimeStepFactor, 1e6))))

    m1 = ArbFloat(0)
    m2 = ArbFloat(0)

    dp  = dpfunc(p, b, dist, m1, m2, G, c, ArbFloat)

    qx1  = - r1
    qy1  = - b/ArbFloat(2)
    qz1  = ArbFloat(0)

    qx2  = r2
    qy2  = b/ArbFloat(2)
    qz2  = ArbFloat(0)

    px1  = p
    py1  = ArbFloat(0)
    pz1  = ArbFloat(0)

    px2  = -p
    py2  = ArbFloat(0)
    pz2  = ArbFloat(0)

    return ( dp , Parameters((false,ArbFloat(1)),
                             (true,ArbFloat("0.005")),
                             (false,ArbFloat("1e-6")),
                             (ArbFloat(0),StartTime+Dur),MaxTimeStep) ,
                    ArbFloat[0;0] ,
                    [qx1;qy1;qz1;qx2;qy2;qz2;px1;py1;pz1;px2;py2;pz2] )
end #-------------------------------------------------------------------

# Test a range of nn values
for nn in [80, 100, 150, 200, 300, 500]
    println("\nTesting with nn = ", nn)
    
    m1  = ArbFloat(0);  # Massless particles
    μ   = ArbFloat(π)/4;
    db  = ArbFloat(100);
    p   = ArbFloat(1);  # Increased momentum for better numerical stability
    b   = ArbFloat(10)*ArbFloat("1.1108305558745590335689712446765042841434478759765625")^nn
    dst   = db*b;
    
    println("Impact parameter b = ", b)
    println("Initial separation dst = ", dst)
    
    mxicd=MXICgenML(p,b,dst;
                        G=ArbFloat(1),c=ArbFloat(1),
                        StartTime=ArbFloat(0),MaxTimeStepFactor=ArbFloat("1e5"))  # Adjusted timestep factor
    
    # Wrap FHE in a form compatible with ODEProblem
    function wrapped_fhe!(du, u, p, t)
        du .= FHE(3, mxicd[3], u)
    end
    
    # Initial timestep and max iterations
    initial_dt = ArbFloat("1e-4")  # Even more conservative timestep for better precision
    max_iterations = Int(1e6)    # Large enough for most simulations
    
    # Define the Hamiltonian gradient function
    function hamiltonian_grad(x)
        return dH(3, mxicd[3], x)
    end

    # Print initial state
    println("Initial state: ", mxicd[4])
    
    # Use hrkintegrator with tcour for adaptive timestepping
    S = hrkintegrator(3, length(mxicd[3]), mxicd[4], 
                      hamiltonian_grad,
                      initial_dt,
                      tcour,
                      mxicd[2].tspan,
                      max_iterations)
    
    # Print final state
    println("Final state: ", S.z[end])
    
    println("Analytical dp (mxicd[1]): ", mxicd[1])
    numerical_dp = S.z[end][11] - S.z[1][11]  # Access first and last states through the z field
    println("Numerical dp: ", numerical_dp)
    println("Percentage difference: ", 100*(abs(numerical_dp)-mxicd[1])/mxicd[1], "%")
end

# hrkintegrator(3, 2, mxicd[4] , x -> pomin.dH_plus_MW(3, mxicd[3], x), δ, no_adapt, mxicd[2].tspan, mxicd[2].iter)
