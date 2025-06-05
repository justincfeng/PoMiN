using LinearAlgebra
using ForwardDiff
using OrdinaryDiffEq
using Logging
using DoubleFloats
using Statistics

include("../pomin.jl")
using .pomin

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
function dpfunc(p, b, d, m1, m2, G, c)
    E1  = c*√((c*m1)^2 + p^2)
    E2  = c*√((c*m2)^2 + p^2)
    dp  = ( (2*G*(E1*E2)^2)/(b*p*(E1+E2)) )*
          ( 1 + ( 1/E1^2 + 1/E2^2 + 4/(E1*E2) )*p^2 + p^4/(E1*E2)^2 )
    return dp
end

"""
    MasslessICgen(p, b, d; G=1, c=1, StartTime=0, MaxTimeStepFactor=1e5)

Generate initial conditions for massless particle collision with:
- p: magnitude of momentum
- b: impact parameter
- d: initial separation
- G: gravitational constant (default 1)
- c: speed of light (default 1)
- StartTime: initial time
- MaxTimeStepFactor: maximum timestep factor

Returns:
- expected_momentum: expected final momentum
- problem_params: ODE problem parameters
- masses: particle masses
- initial_state: initial phase space state
"""
function MasslessICgen(p, b, d;
    G=Double64(1), c=Double64(1),
    StartTime=Double64(0), MaxTimeStepFactor=Double64(1e5))

    # Both particles are massless
    m1 = Double64(0.0)
    m2 = Double64(0.0)
    masses = [m1, m2]

    # Initial positions
    # Particle 1 starts at (-d/2, -b/2, 0)
    # Particle 2 starts at (d/2, b/2, 0)
    qx1 = Double64(-d/2)
    qy1 = Double64(-b/2)
    qz1 = Double64(0)
    
    qx2 = Double64(d/2)
    qy2 = Double64(b/2)
    qz2 = Double64(0)

    # Initial momenta
    # Particle 1 moves right with momentum p
    # Particle 2 moves left with momentum p
    px1 = Double64(p)
    py1 = Double64(0)
    pz1 = Double64(0)
    
    px2 = Double64(-p)
    py2 = Double64(0)
    pz2 = Double64(0)

    # Initial state vector
    Z_init = [qx1;qy1;qz1;qx2;qy2;qz2;px1;py1;pz1;px2;py2;pz2]

    # Expected momentum conservation
    expected_momentum = p

    # Integration timespan
    tmax = MaxTimeStepFactor * d/c
    tspan = (StartTime, tmax)

    # Problem parameters
    problem_params = (
        tspan = tspan,
        expected_momentum = expected_momentum
    )

    return expected_momentum, problem_params, masses, Z_init
end

# Test parameters
nn = 250  # Number of iterations
p = Double64(1.0)  # Initial momentum
db = Double64(100000)  # Base distance factor
b_base = Double64(1.1108305558745590335689712446765042841434478759765625)

# Run tests with increasing impact parameter
for i in 1:nn
    b = Double64(10) * b_base^i  # Impact parameter grows exponentially
    d = db * b  # Initial separation scales with impact parameter
    
    # Generate initial conditions
    mdata = MasslessICgen(p, b, d;
        G=Double64(1), c=Double64(1),
        StartTime=Double64(0), MaxTimeStepFactor=Double64(1e5))
    
    # Integrate using high precision
    sol = pomin.jlintegratorfull(
        mdata[4],  # Initial state
        (du,u,p,t) -> pomin.Jsympl(pomin.dH(3, mdata[3], u)),  # EOM
        mdata[2].tspan,  # Time span
        mdata[3],  # Masses
        Double64(1e-25),  # abstol
        Double64(1e-25),  # reltol
        AutoTsit5(Rodas5())  # Integrator
    )
    
    # Extract initial and final momenta
    initial_p1 = mdata[4][7:9]   # Initial p1 vector
    initial_p2 = mdata[4][10:12] # Initial p2 vector
    final_p1 = sol[end][7:9]     # Final p1 vector
    final_p2 = sol[end][10:12]   # Final p2 vector
    
    # Calculate scattering angles
    angle1 = acos(dot(initial_p1, final_p1)/(norm(initial_p1)*norm(final_p1)))
    angle2 = acos(dot(initial_p2, final_p2)/(norm(initial_p2)*norm(final_p2)))
    
    # Calculate transverse momentum changes
    dp_num1 = norm(initial_p1) * sin(angle1)
    dp_num2 = norm(initial_p2) * sin(angle2)
    
    # Calculate analytical momentum change
    dp_analytical = dpfunc(p, b, d, Double64(0.0), Double64(0.0), Double64(1.0), Double64(1.0))
    
    # Calculate relative errors
    error1 = abs(dp_num1 - dp_analytical)/dp_analytical
    error2 = abs(dp_num2 - dp_analytical)/dp_analytical
    
    # Debug output
    println("  Initial p1: $(initial_p1)")
    println("  Initial p2: $(initial_p2)")
    println("  Final p1: $(final_p1)")
    println("  Final p2: $(final_p2)")
    println("  Scattering angle 1: $(angle1*180/π)°")
    println("  Scattering angle 2: $(angle2*180/π)°")
    
    println("Test $i: b = $b")
    println("  Initial separation: $(d)")
    println("  Analytical dp: $(dp_analytical)")
    println("  Numerical dp1: $(dp_num1)")
    println("  Numerical dp2: $(dp_num2)")
    println("  Relative error 1: $(100*error1)%")
    println("  Relative error 2: $(100*error2)%")
    println("  Integration time: $(mdata[2].tspan[2])")
    println()
end
