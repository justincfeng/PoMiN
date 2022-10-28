using DoubleFloats

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

function MXICgen(m1,μ,p,b,d;
                    G=Double64(1),c=Double64(1),
                    StartTime=Double64(0),MaxTimeStepFactor=Double64(1e5))

    n   = Double64(2)

    if m1 == Double64(0.)
        Dur = d/c
        r1  = d/2
        r2  = d/2
    else
        v1  = p*c/√((c*m1)^2 + p^2)
        if μ == Double64(0.)
            τ   = d/(v1+c)
            Dur = 2*τ
            r1  = v1*τ
            r2  = c*τ
        else
            v2  = p*c/√((μ*c*m1)^2 + p^2)
            τ   = d/(v1+v2)
            Dur = 2*τ
            r1  = v1*τ
            r2  = v2*τ  
        end
    end

    MaxTimeStep = Int(Dur/MaxTimeStepFactor)

    E1  = c*√((c*m1)^2 + p^2)
    E2  = c*√((μ*c*m1)^2 + p^2)

    dp  = (2*G*(E1*E2)^2)*
          (
            Double64(1.)
            + ( 1/E1^2 + 1/E2^2 + 4/(E1*E2) )*p^2
            + p^4/(E1*E2)^2
          )
          /
          (b*p*(E1+E2))

    dpm = 4*G*p^2 / b

    m2  = μ*m1

    qx1  = - r1
    qy1  = - b/2
    qz1  = Double64(0.)

    qx2  = r2
    qy2  = b/2
    qz2  = Double64(0.)

    px1  = p
    py1  = Double64(0.)
    pz1  = Double64(0.)

    px2  = p
    py2  = Double64(0.)
    pz2  = Double64(0.)

    return (
    dp , 
    Parameters((false,1.),(true,0.005),(false,1e-6),(0.,StartTime+Dur),MaxTimeStep) ,
    [qx1;qy1;qz1;qx2;qy2;qz2;px1;py1;pz1;px2;py2;pz2]
    )
end

m1  = 0.0497992;
μ   = π/4;
db  = 100000;
p   = 10*m1;
b   = 1.110830555874559*10.;
d   = db*b;

mxicd=MXICgen(m1,μ,p,b,d;
                    G=Double64(1),c=Double64(1),
                    StartTime=Double64(0),MaxTimeStepFactor=Double64(1e5))
