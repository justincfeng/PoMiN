using DoubleFloats

include("../broyden.jl")

function closestApproachStraightLine(v_init::RealVec)

    masses = Double64[1]        # one object, the Sun, with mass in units of solar masses
    tspan = (Double64(0), Double64(3E+13))
    δ = Double64(3E+10)
    no_adapt = (δ, z, zdot) -> δ  # turns off adaptive time stepping
    maxit = 1000
    p_init = masses[1] * v_init
    Z_init = [0, 0, 0, p_init[1], p_init[2], p_init[3]]
    sol = hrkintegrator(3, 1, Z_init, x -> pomin.dH(3, masses, x), δ, no_adapt, tspan, maxit)

    # get final position of particle 1
    numsteps = length(sol.t)
    q = Z2q(sol.N, sol.d, 1, sol.z[numsteps])

    return q - Double64[21213203435596.42729, 0, 0]
end

function closestApproach2masses(v_init::RealVec)

    masses = Double64[1, 1]
    tspan = (Double64(0), Double64(3E+13))
    δ = Double64(3E+10)
    no_adapt = (δ, z, zdot) -> δ  # turns off adaptive time stepping
    maxit = 1000
    p_init = masses[1] * v_init
    Z_init = [0, 0, 0, p_init[1], p_init[2], p_init[3], 0, 1e8, 0, 0, 0, 0]
    sol = hrkintegrator(3, 1, Z_init, x -> pomin.dH(3, masses, x), δ, no_adapt, tspan, maxit)

    # get final position of particle 1
    numsteps = length(sol.t)
    q = Z2q(sol.N, sol.d, 1, sol.z[numsteps])

    return q - Double64[-1.1997161984344225e14, -2.9113501694797115e9, 0.0]
end

function BroydenMethod(v_0)
    J = ForwardDiff.jacobian(x -> closestApproach2masses(x), v_0)
    println("J=",J)
    f_0 = closestApproach2masses(v_0)
    println(stderr, "Init norm F = ", norm(f_0))
    return bsolve(x -> closestApproach2masses(x), J, f_0, v_0)
end

v_init = [0.5, 0.5, 0]
println(BroydenMethod(v_init))

# println(closestApproach2masses([0.8,0,0]))