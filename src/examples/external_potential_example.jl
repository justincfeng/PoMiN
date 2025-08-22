

using DoubleFloats
using LinearAlgebra
using Printf

# Include the necessary PoMiN modules
include("../pomin.jl")
include("../core/pomin-types.jl")  # Load type definitions first
include("../core/physics/Hamiltonians/HamTools.jl")  # For Z2q function
include("../core/physics/external_potentials/external.jl")

using Plots

# Import the gradient operator from pomin module
∂ = pomin.∂

"""
    elliptical_initial_conditions(a,m,e,periapsis=false)
"""
function elliptical_initial_conditions(a,M,e,periapsis=false)
    if periapsis
        r0 = a*(1-e)  # Periapsis
    else
        r0 = a*(1+e)  # Apoapsis
    end
    ve = sqrt(abs( M*(2.0/r0-1.0/a) ))
    l  = ve*r0            # Angular momentum
    ϵ = ve^2/2 - M/r0     # Total energy
    Te = 2*π*sqrt(a^3/M)
    return (r0,ve,l,ϵ,Te)
end #-------------------------------------------------------------------

r0,ve,l,ϵ,Te = elliptical_initial_conditions(10.0,1.0,0.1)

V = ΦCentralMass( zeros(3) , Double64 , 1.0 )

m = [1e-10,1e-10]

U = UConstructor(V,m,zeros(3))

q0 = [[r0,0.0,0.0],[-r0,0.0,0.0]]  # Position vector as array of vectors
p0 = [[0.0,m[1]*ve,0.0],[0.0,-m[2]*ve,0.0]]  # Momentum vector as array of vectors

system = pomin.Particles(m, q0, p0)

tspan = (0.0, Te)  # From t=0 to one orbital period
params = pomin.ParametersJulia(tspan)

sol = pomin.solve(system, params, U)

#   PLOTTING

# Extract trajectory data for plotting
n_particles = length(m)
d = 3  # 3D space
n_steps = length(sol.t)

# Extract positions for each particle over time
x1_traj = [sol.u[i][1] for i in 1:n_steps]  # Particle 1 x-position
y1_traj = [sol.u[i][2] for i in 1:n_steps]  # Particle 1 y-position
z1_traj = [sol.u[i][3] for i in 1:n_steps]  # Particle 1 z-position

x2_traj = [sol.u[i][4] for i in 1:n_steps]  # Particle 2 x-position
y2_traj = [sol.u[i][5] for i in 1:n_steps]  # Particle 2 y-position
z2_traj = [sol.u[i][6] for i in 1:n_steps]  # Particle 2 z-position

# Create orbital trajectory plot
println("Creating orbital trajectory plots...")

# 2D projection in x-y plane
p1 = plot(x1_traj, y1_traj, 
          label="Particle 1", 
          linewidth=2, 
          color=:blue,
          title="Binary System in External Central Mass Potential",
          xlabel="x", ylabel="y",
          aspect_ratio=:equal)
plot!(p1, x2_traj, y2_traj, 
      label="Particle 2", 
      linewidth=2, 
      color=:red)

# Mark initial positions
scatter!(p1, [x1_traj[1]], [y1_traj[1]], 
         label="Start 1", 
         markersize=6, 
         color=:blue, 
         markershape=:circle)
scatter!(p1, [x2_traj[1]], [y2_traj[1]], 
         label="Start 2", 
         markersize=6, 
         color=:red, 
         markershape=:circle)

# Mark central mass at origin
scatter!(p1, [0], [0], 
         label="Central Mass", 
         markersize=8, 
         color=:black, 
         markershape=:star)

# 3D trajectory plot
p2 = plot(x1_traj, y1_traj, z1_traj,
          label="Particle 1",
          linewidth=2,
          color=:blue,
          title="3D Binary Orbits",
          xlabel="x", ylabel="y", zlabel="z")
plot!(p2, x2_traj, y2_traj, z2_traj,
      label="Particle 2",
      linewidth=2,
      color=:red)

# Mark central mass at origin in 3D
scatter!(p2, [0], [0], [0],
         label="Central Mass",
         markersize=8,
         color=:black,
         markershape=:star)

# Time evolution of separation
separation = [sqrt((x1_traj[i] - x2_traj[i])^2 + 
                  (y1_traj[i] - y2_traj[i])^2 + 
                  (z1_traj[i] - z2_traj[i])^2) for i in 1:n_steps]

p3 = plot(sol.t, separation,
          label="Binary Separation",
          linewidth=2,
          color=:green,
          title="Binary Separation vs Time",
          xlabel="Time", ylabel="Separation")

# Combine plots
combined_plot = plot(p1, p2, p3, layout=(2,2), size=(800,600))

# Display the plot
display(combined_plot)

# Save the plot
savefig(combined_plot, "binary_orbits_external_potential.png")
println("Plot saved as 'binary_orbits_external_potential.png'")

# Print some orbital statistics
println("\n=== Orbital Statistics ===")
println("Integration time: $(sol.t[end]) time units")
println("Number of time steps: $n_steps")
println("Initial separation: $(separation[1])")
println("Final separation: $(separation[end])")
println("Max separation: $(maximum(separation))")
println("Min separation: $(minimum(separation))")
println("Average separation: $(sum(separation)/length(separation))")
