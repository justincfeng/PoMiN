#!/usr/bin/env julia

# Simple scattering example using PoMiN

# Add the parent directory to the load path to find the module
pushfirst!(LOAD_PATH, dirname(@__DIR__))
pushfirst!(LOAD_PATH, joinpath(dirname(@__DIR__), "src"))

# Import required packages
using LinearAlgebra

# Import the pomin module
include(joinpath(dirname(@__DIR__), "src", "pomin.jl"))
using .pomin

# Parameters for the scattering simulation
# -----------------------------------------
# Masses of the two particles
m1 = 1.0
m2 = 1.0

# Initial momentum in center of mass frame
p_init = 0.1

# Impact parameter (perpendicular distance between trajectories)
impact_param = 10.0

# Initial separation distance along x-axis
init_separation = 50.0

# Speed of light (in natural units)
c = 1.0

# Simulation parameters
# -----------------------------------------
# Time span: from t_start to t_end
t_start = 0.0
t_end = 200.0

# Create the particle system for scattering
system = setup_scattering(m1, m2, p_init, impact_param, init_separation, c=c)

# Print initial configuration
println("Initial configuration:")
println("Particle 1: mass = $(system.m[1]), position = $(system.q[1]), momentum = $(system.p[1])")
println("Particle 2: mass = $(system.m[2]), position = $(system.q[2]), momentum = $(system.p[2])")

# Create simulation parameters
# Use the Julia integrator with Tsit5 algorithm
params = Parameters(
    # rkl: Use RK4 integrator? (false, Courant number)
    (false, 0.1),
    # jli: Use Julia integrator? (true, tolerance, integrator method)
    (true, 1e-10, "Tsit5"),
    # tspan: Time span (start, end)
    (t_start, t_end),
    # iter: Number of iterations
    1000
)

# Run the simulation
println("\nRunning simulation...")
solution = solve(system, params)

# Extract results
# The solution contains trajectories for all particles
times = solution.t
positions = []
trajectories = []

# Process the solution to extract particle trajectories
for i in 1:length(times)
    # Each solution point contains the full state vector
    state = solution.u[i]
    
    # Extract positions for each particle (assuming 3D)
    dim = length(system.q[1])
    num_particles = length(system.m)
    
    # For this time step, get positions of all particles
    pos_at_time = []
    for j in 1:num_particles
        # Extract position of particle j (first dim*num_particles elements are positions)
        pos_j = state[(j-1)*dim+1:j*dim]
        push!(pos_at_time, pos_j)
    end
    
    push!(positions, pos_at_time)
    
    # Store as trajectories for plotting
    if i == 1
        for j in 1:num_particles
            push!(trajectories, [positions[i][j]])
        end
    else
        for j in 1:num_particles
            push!(trajectories[j], positions[i][j])
        end
    end
end

# Print final configuration
println("\nFinal configuration:")
final_pos = positions[end]
final_state = solution.u[end]
dim = length(system.q[1])
num_particles = length(system.m)

for j in 1:num_particles
    pos_j = final_state[(j-1)*dim+1:j*dim]
    mom_j = final_state[num_particles*dim + (j-1)*dim+1:num_particles*dim + j*dim]
    println("Particle $j: position = $pos_j, momentum = $mom_j")
end

# Calculate deflection angles
println("\nCalculating deflection angles...")

# Extract initial and final momenta
initial_p1 = system.p[1]
initial_p2 = system.p[2]

final_state = solution.u[end]
dim = length(system.q[1])
num_particles = length(system.m)

final_p1 = final_state[num_particles*dim + 1:num_particles*dim + dim]
final_p2 = final_state[num_particles*dim + dim+1:num_particles*dim + 2*dim]

# Calculate deflection angles
deflection1 = acos(dot(initial_p1, final_p1)/(norm(initial_p1)*norm(final_p1))) * 180/π
deflection2 = acos(dot(initial_p2, final_p2)/(norm(initial_p2)*norm(final_p2))) * 180/π

println("Deflection angle for particle 1: $(deflection1) degrees")
println("Deflection angle for particle 2: $(deflection2) degrees")

println("\nSimulation completed successfully!")
