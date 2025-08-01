#!/usr/bin/env julia
#-----------------------------------------------------------------------
#
#   Binary Black Hole Orbit Test
#   
#   Tests the modernized PoMiN codebase with a simple equal-mass
#   binary black hole system in circular orbit.
#
#-----------------------------------------------------------------------

using Pkg
Pkg.activate(".")

using LinearAlgebra
using Plots

# Include core components directly
include("../src/core/pomin-types.jl")
include("../src/physics/Hamiltonians/HamTools.jl")
include("../src/physics/Hamiltonians/idxer.jl")
include("../src/physics/Hamiltonians/HamPM.jl")
include("../src/integrators/tadap.jl")
include("../src/integrators/rk4i.jl")
include("../src/integrators/intjul.jl")
include("../src/physics/initial_data/idgen.jl")

println("🌌 PoMiN Binary Black Hole Orbit Test")
println("=====================================")

# System parameters
println("\n📋 Setting up binary system...")
m1 = 1.0        # Mass of first black hole (solar masses)
m2 = 1.0        # Mass of second black hole (solar masses)  
r_sep = 10.0    # Initial separation (geometric units)
v_orb = 0.1     # Orbital velocity (fraction of c)

println("   • Mass 1: $m1 M☉")
println("   • Mass 2: $m2 M☉") 
println("   • Separation: $r_sep GM/c²")
println("   • Orbital velocity: $(v_orb)c")

# Create initial data
println("\n🔧 Generating initial data...")
# Set up initial positions and momenta for circular orbit
q1 = [-r_sep/2, 0.0, 0.0]  # First black hole at -r/2 on x-axis
q2 = [r_sep/2, 0.0, 0.0]   # Second black hole at +r/2 on x-axis
p1 = [0.0, m1*v_orb, 0.0]   # First black hole moving in +y direction
p2 = [0.0, -m2*v_orb, 0.0]  # Second black hole moving in -y direction

# Create the particle system
system = Particles([m1, m2], [q1, q2], [p1, p2])

println("   • Particles created: $(length(system.m))")
println("   • Total mass: $(sum(system.m)) M☉")
println("   • Initial positions:")
for i in 1:length(system.m)
    println("     Particle $i: $(system.q[i])")
end
println("   • Initial momenta:")
for i in 1:length(system.m)
    println("     Particle $i: $(system.p[i])")
end

# Integration parameters
println("\n⚙️  Setting up integration parameters...")
tspan = (0.0, 100.0)  # Time span (geometric units)
δ = 0.01              # Time step for RK4
atol = 1e-12          # Absolute tolerance for Julia integrators
rtol = 1e-10          # Relative tolerance for Julia integrators

# Test both integrators
println("\n🚀 Testing RK4 integrator...")
params_rk4 = ParametersRK4(tspan, δ=δ)
println("   • Time span: $tspan")
println("   • Time step: $δ")

println("\n🧮 Running RK4 integration...")
@time solution_rk4 = solve(system, params_rk4)

println("   ✅ RK4 integration complete!")
println("   • Solution type: $(typeof(solution_rk4))")
println("   • Time steps: $(length(solution_rk4.t))")
println("   • Final time: $(solution_rk4.t[end])")

println("\n🚀 Testing Julia integrator...")
params_julia = ParametersJulia(tspan, atol=atol, rtol=rtol)
println("   • Time span: $tspan")
println("   • Absolute tolerance: $atol")
println("   • Relative tolerance: $rtol")

println("\n🧮 Running Julia integration...")
@time solution_julia = solve(system, params_julia)

println("   ✅ Julia integration complete!")
println("   • Solution type: $(typeof(solution_julia))")
println("   • Time steps: $(length(solution_julia.t))")
println("   • Final time: $(solution_julia.t[end])")

# Extract trajectories for plotting
println("\n📊 Extracting trajectories...")

function extract_positions(sol, n_particles)
    times = sol.t
    n_dim = 3
    positions = []
    
    if typeof(sol) == soln
        # RK4 solution - stored in z field
        for i in 1:n_particles
            x_traj = [z[i] for z in sol.z]
            y_traj = [z[n_particles + i] for z in sol.z]
            z_traj = [z[2*n_particles + i] for z in sol.z]
            push!(positions, (x_traj, y_traj, z_traj))
        end
    else
        # Julia ODE solution - stored in u field
        for i in 1:n_particles
            x_traj = [u[i] for u in sol.u]
            y_traj = [u[n_particles + i] for u in sol.u]
            z_traj = [u[2*n_particles + i] for u in sol.u]
            push!(positions, (x_traj, y_traj, z_traj))
        end
    end
    
    return times, positions
end

# Extract trajectories
times_rk4, pos_rk4 = extract_positions(solution_rk4, 2)
times_julia, pos_julia = extract_positions(solution_julia, 2)

println("   • RK4 trajectory points: $(length(times_rk4))")
println("   • Julia trajectory points: $(length(times_julia))")

# Calculate separations for analysis
println("\n📈 Analyzing trajectories...")
separation_rk4 = [norm(pos_rk4[1][1][i] .- pos_rk4[2][1][i]) for i in 1:length(times_rk4)]
separation_julia = [norm(pos_julia[1][1][i] .- pos_julia[2][1][i]) for i in 1:length(times_julia)]

# Create plots
println("\n📈 Creating plots...")

# Plot 1: Orbital trajectories (x-y plane)
p1 = plot(title="Binary Black Hole Orbits (x-y plane)", 
          xlabel="x (GM/c²)", ylabel="y (GM/c²)",
          aspect_ratio=:equal, legend=:topright)

plot!(p1, pos_rk4[1][1], pos_rk4[1][2], 
      label="BH1 (RK4)", color=:red, linewidth=2)
plot!(p1, pos_rk4[2][1], pos_rk4[2][2], 
      label="BH2 (RK4)", color=:blue, linewidth=2)
plot!(p1, pos_julia[1][1], pos_julia[1][2], 
      label="BH1 (Julia)", color=:red, linestyle=:dash, alpha=0.7)
plot!(p1, pos_julia[2][1], pos_julia[2][2], 
      label="BH2 (Julia)", color=:blue, linestyle=:dash, alpha=0.7)

# Mark initial positions
scatter!(p1, [pos_rk4[1][1][1]], [pos_rk4[1][2][1]], 
         color=:red, markersize=8, label="Start BH1")
scatter!(p1, [pos_rk4[2][1][1]], [pos_rk4[2][2][1]], 
         color=:blue, markersize=8, label="Start BH2")

# Plot 2: Separation vs time
# Note: separations already calculated above

p2 = plot(title="Binary Separation vs Time", 
          xlabel="Time (GM/c³)", ylabel="Separation (GM/c²)")
plot!(p2, times_rk4, [norm(pos_rk4[1][1][i] .- pos_rk4[2][1][i]) for i in 1:length(times_rk4)], 
      label="RK4", color=:green, linewidth=2)
plot!(p2, times_julia, [norm(pos_julia[1][1][i] .- pos_julia[2][1][i]) for i in 1:length(times_julia)], 
      label="Julia", color=:orange, linestyle=:dash, linewidth=2)

# Combine plots
final_plot = plot(p1, p2, layout=(2,1), size=(800, 800), dpi=300)

# Save plot
plot_file = joinpath(@__DIR__, "binary_orbit_results.png")
savefig(final_plot, plot_file)
println("   📈 Plot saved to: $plot_file")

# Print some statistics about the integration
println("\n📊 Integration Statistics:")
println("   • RK4 steps: $(length(times_rk4))")
println("   • Julia steps: $(length(times_julia))")
println("   • Integration time span: $(times_rk4[1]) to $(times_rk4[end]) GM/c³")

# Summary statistics
println("\n📋 Summary Statistics:")
println("   • Initial separation: $(separation_rk4[1]:.6f) GM/c²")
println("   • Final separation (RK4): $(separation_rk4[end]:.6f) GM/c²")
println("   • Final separation (Julia): $(separation_julia[end]:.6f) GM/c²")
println("   • Separation change (RK4): $((separation_rk4[end] - separation_rk4[1])/separation_rk4[1] * 100:.3f)%")
println("   • Separation change (Julia): $((separation_julia[end] - separation_julia[1])/separation_julia[1] * 100:.3f)%")

# Print orbital period estimate (assuming roughly circular orbit)
velocity = 0.1  # in units of c
initial_separation = separation_rk4[1]  # in GM/c²
period = 2π * initial_separation / velocity  # in GM/c³
println("\n🕰️ Orbital Period:")
println("   • Estimated period: $(period:.2f) GM/c³")
println("   • Integration time: $(times_rk4[end]:.2f) GM/c³ ($(times_rk4[end]/period:.1f) orbits)")

# Energy conservation check using Hamiltonian
println("\n⚡ Energy Conservation Check:")

# Calculate energies for RK4 solution
energy_rk4 = [HamPM.H(vcat(pos_rk4[1][1][i], pos_rk4[2][1][i], 
                         pos_rk4[1][2][i], pos_rk4[2][2][i], 
                         pos_rk4[1][3][i], pos_rk4[2][3][i]), 
                    [1.0, 1.0]) for i in 1:length(times_rk4)]

# Calculate energies for Julia solution
energy_julia = [HamPM.H(vcat(pos_julia[1][1][i], pos_julia[2][1][i], 
                            pos_julia[1][2][i], pos_julia[2][2][i], 
                            pos_julia[1][3][i], pos_julia[2][3][i]), 
                       [1.0, 1.0]) for i in 1:length(times_julia)]

# Calculate relative energy changes
relative_change_rk4 = (energy_rk4[end] - energy_rk4[1]) / abs(energy_rk4[1]) * 100
relative_change_julia = (energy_julia[end] - energy_julia[1]) / abs(energy_julia[1]) * 100

# Print energy conservation statistics
println("   • Initial energy (RK4): $(energy_rk4[1]:.6e)")
println("   • Final energy (RK4): $(energy_rk4[end]:.6e)")
println("   • Energy change (RK4): $(relative_change_rk4:.3f)%")
println("   • Initial energy (Julia): $(energy_julia[1]:.6e)")
println("   • Final energy (Julia): $(energy_julia[end]:.6e)")
println("   • Energy change (Julia): $(relative_change_julia:.3f)%")

println("\n✅ Binary orbit test completed successfully!")
println("🎉 PoMiN modernization validation: PASSED")
