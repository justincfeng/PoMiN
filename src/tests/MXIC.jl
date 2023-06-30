using DoubleFloats
using LinearAlgebra
using ForwardDiff
using OrdinaryDiffEq

include("../pomin-types.jl")        # data type definitions
include("../pomin-io.jl")           # input/output routines
include("../idxer.jl")              # routines for index management
include("../integrators/rk4i.jl")   # Julia integrator functions
include("../HamPM.jl")              # Julia integrator functions
include("../integrators/int.jl")    # Julia integrator functions

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

#-----------------------------------------------------------------------
function MXICgen(m1,μ,p,b,d;
                 G=Double64(1),c=Double64(1),
                 StartTime=Double64(0),MaxTimeStepFactor=Double64(1e5))

    n   = Double64(2)

    if m1 == Double64(0.)
        Dur = d/c
        r1  = d/Double64(2)
        r2  = d/Double64(2)
    else
        v1  = p*c/√((c*m1)^2 + p^2)
        if μ == Double64(0.)
            τ   = d/(v1+c)
            Dur = Double64(2)*τ
            r1  = v1*τ
            r2  = c*τ
        else
            v2  = p*c/√((μ*c*m1)^2 + p^2)
            τ   = d/(v1+v2)
            Dur = Double64(2)*τ
            r1  = v1*τ
            r2  = v2*τ  
        end
    end

    Dur = 2*d*(p/m1)

    MaxTimeStep = Int(Dur/MaxTimeStepFactor)

    E1  = c*√((c*m1)^2 + p^2)
    E2  = c*√((μ*c*m1)^2 + p^2)

    dp  = ( (Double64(2)*G*(E1*E2)^2)/(b*p*(E1+E2)) )*
          ( Double64(1) + 
            ( Double64(1)/E1^2 + Double64(1)/E2^2 + 
              Double64(4)/(E1*E2) )*p^2 + p^4/(E1*E2)^2 )

    #dpm = 4*G*p^2 / b

    m2  = μ*m1

    #qx1  = - r1
    qx1  = -d
    qy1  = - b/Double64(2)
    qz1  = Double64(0.)

    #qx2  = r2
    qx2  = d
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
                    [m1;m2] ,
                    [qx1;qy1;qz1;qx2;qy2;qz2;px1;py1;pz1;px2;py2;pz2] )
end #-------------------------------------------------------------------

m1  = Double64(0.0497992);
μ   = Double64(π/4);
db  = Double64(100000);
p   = Double64(10*m1);
b   = Double64(10*(1.1108305558745590335689712446765042841434478759765625)^2);
d   = Double64(db*b);

mxicd=MXICgen(m1,μ,p,b,d;
                    G=Double64(1),c=Double64(1),
                    StartTime=Double64(0),MaxTimeStepFactor=Double64(1e5))


S = hointegrator( mxicd[4] , z->FHE(3,mxicd[3],z) , mxicd[2].tspan )

# hrkintegrator(3, 2, mxicd[4] , x -> pomin.dH_plus_MW(3, mxicd[3], x), δ, no_adapt, mxicd[2].tspan, mxicd[2].iter)

#-----------------------------------------------------------------------
function MXICmassive(m1,μ,p,b,d;
    G=Double64(1),c=Double64(1),
    StartTime=Double64(0),MaxTimeStepFactor=Double64(1e5))

    hrkintegrator(3, 2, Z_init, x -> pomin.dH_plus_MW(3, masses, x), δ, no_adapt, tspan, maxit)

end #-------------------------------------------------------------------