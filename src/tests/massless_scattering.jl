using DoubleFloats
using LinearAlgebra
using ForwardDiff
using OrdinaryDiffEq

include("../pomin-types.jl")        # data type definitions
include("../pomin-io.jl")           # input/output routines
include("../idxer.jl")              # routines for index management
include("../integrators/Tao.jl")    # Symplectic integrator
include("../HamPM.jl")              # Hamiltonian functions
include("../integrators/intjul.jl") # Julia integrator functions

#   Datatype for parameters
struct Parameters
    sym::Tuple{Bool,Real}     # Use the symplectic integrator?            The real parameter in the tuple is the ω parameter in the Tao map
    rkl::Tuple{Bool,Real}     # Use the rk4 integrator?                   The real parameter in the tuple is the Courant number
    jli::Tuple{Bool,Real}     # Use integrators in OrdinaryDiffEq.jl?     The real parameter in the tuple is the tolerance
    tspan::Tuple{Real,Real}   # Time span (initial_time,final_time)
    iter::Int                # Either the number of iterations, or maximum number of iterations, depending on integrator
end

# Both particles are massless
m1 = Double64(0.0)
m2 = Double64(0.0)

# Initial momentum magnitude (same for both particles)
p = Double64(1/2)

# Impact parameter and initial separation
b = Double64(6e5)        # Impact parameter (reduced further)
d = Double64(b*10^2)      # Initial separation (reduced further)

#-----------------------------------------------------------------------
function MasslessScatteringGen(p, b, d;
                              G=Float64(1), c=Double64(1),
                              StartTime=Double64(0), MaxTimeStepFactor=Double64(1e8))  # Increased time resolution further
    
    # For massless particles, velocity = c
    Dur = 2*d/c  # Time for particles to travel total distance 2d
    
    # Energy for massless particles is just |p|
    E1 = c*p
    E2 = c*p
    
    # Scattering angle estimate from post-Minkowskian approximation
    dp = ((Double64(2)*G*(E1*E2)^2)/(b*p*(E1+E2))) * 
         (Double64(1) + Double64(4)/(E1*E2) + p^2/(E1*E2)^2)
    
    MaxTimeStep = Int(Dur/MaxTimeStepFactor)
    
    # Initial positions - symmetric setup
    qx1 = -d
    qy1 = -b/Double64(2)
    qz1 = Double64(0.)
    
    qx2 = d
    qy2 = b/Double64(2)
    qz2 = Double64(0.)
    
    # Initial momenta - equal and opposite
    px1 = p
    py1 = Double64(0.)
    pz1 = Double64(0.)
    
    px2 = -p
    py2 = Double64(0.)
    pz2 = Double64(0.)
    
    return (dp, Parameters((false, Double64(1.)),
                          (true, Double64(0.00001)),  # Reduced time step further
                          (false, Double64(1e-6)),
                          (Double64(0.), StartTime+Dur), MaxTimeStep),
            [m1; m2],
            [qx1; qy1; qz1; qx2; qy2; qz2; px1; py1; pz1; px2; py2; pz2])
end #-------------------------------------------------------------------

# Number of test points
nn = 250

# Generate initial conditions
dp, par, mass, Z = MasslessScatteringGen(p, b, d)

# Run simulation
function F!(du, u, p, t)
    du .= dH_plus_MW(3, mass, u)
end

# Calculate and print conserved quantities
function calculate_energy(z)
    p1 = sqrt(z[7]^2 + z[8]^2 + z[9]^2)
    p2 = sqrt(z[10]^2 + z[11]^2 + z[12]^2)
    return p1 + p2  # For massless particles, E = |p|
end

function calculate_momentum(z)
    return [z[7] + z[10], z[8] + z[11], z[9] + z[12]]
end

function compute_hamiltonian(z, masses)
    # Extract positions and momenta
    x1 = z[1:3]
    x2 = z[4:6]
    p1 = z[7:9]
    p2 = z[10:12]
    
    # Compute H2 term (kinetic)
    H2 = sqrt(dot(p1,p1)) + sqrt(dot(p2,p2))
    
    # Compute separation
    r = sqrt(sum((x1 .- x2).^2))
    
    # Compute H3 term (potential)
    H3 = -2.0/(r)
    
    return H2, H3, H2 + H3
end

function print_state(t, z)
    p1 = sqrt(z[7]^2 + z[8]^2 + z[9]^2)
    p2 = sqrt(z[10]^2 + z[11]^2 + z[12]^2)
    r = sqrt((z[1]-z[4])^2 + (z[2]-z[5])^2 + (z[3]-z[6])^2)
    
    # Calculate velocities (for massless particles, v = p/|p|)
    v1 = [z[7]/p1, z[8]/p1, z[9]/p1]
    v2 = [z[10]/p2, z[11]/p2, z[12]/p2]
    rel_v = sqrt(sum((v1 .- v2).^2))
    
    # Calculate Hamiltonian terms
    H2, H3, H_total = compute_hamiltonian(z, [0.0, 0.0])
    
    println("Time: $t")
    println("Separation: $r")
    println("Particle 1: |p| = $p1, v = $v1")
    println("Particle 2: |p| = $p2, v = $v2")
    println("Relative velocity: $rel_v")
    println("Total momentum: $(calculate_momentum(z))")
    println("H2 (kinetic): $H2")
    println("H3 (potential): $H3")
    println("H_total: $H_total")
    println("-----------------")
end

initial_energy = calculate_energy(Z)
initial_momentum = calculate_momentum(Z)

# Print initial state
println("Initial state:")
print_state(par.tspan[1], Z)

# Print initial Hamiltonian
println("\nInitial Hamiltonian analysis:")
H2_init, H3_init, H_total_init = compute_hamiltonian(Z, [0.0, 0.0])
println("Initial H2 (kinetic): $H2_init")
println("Initial H3 (potential): $H3_init")
println("Initial H_total: $H_total_init")
println("-----------------")

# Run simulation with more frequent output
if par.rkl[1]
    prob = ODEProblem(F!, Z, par.tspan)
    sol = solve(prob, RK4(), dt=par.rkl[2], saveat=LinRange(par.tspan[1], par.tspan[2], 50))
    t = sol.t
    z = sol.u
    
    # Print intermediate states
    for i in 1:length(t)
        print_state(t[i], z[i])
    end
elseif par.sym[1]
    t, z = tao(mass, Z, par.tspan[1], par.tspan[2], par.iter, par.sym[2])
elseif par.jli[1]
    t, z = jlintegrator(mass, Z, par.tspan[1], par.tspan[2], par.jli[2])
end

# Print final state
println("\nFinal state:")
print_state(t[end], z[end])

final_energy = calculate_energy(z[end])
final_momentum = calculate_momentum(z[end])

println("\nFinal Hamiltonian analysis:")
H2_final, H3_final, H_total_final = compute_hamiltonian(sol.u[end], [0.0, 0.0])
println("Final H2 (kinetic): $H2_final")
println("Final H3 (potential): $H3_final")
println("Final H_total: $H_total_final")
println("Relative Hamiltonian error: $(abs(H_total_final - H_total_init)/abs(H_total_init))")
println("-----------------")

println("\nConservation analysis:")
println("Energy:")
println("  Initial: $initial_energy")
println("  Final:   $final_energy")
println("  Relative error: $(abs(final_energy - initial_energy)/initial_energy)")
println("\nMomentum:")
println("  Initial: $initial_momentum")
println("  Final:   $final_momentum")
println("  Relative error: $(norm(final_momentum - initial_momentum)/norm(initial_momentum))")
