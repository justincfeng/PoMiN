

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
    milky_way_orbit_parameters(V_potential)
    
Calculates orbital parameters using the actual MW potential gradient at 8.4 kpc.
Uses proper solar mass units where 1 M_sun = 1476.67 m and c = 1.
"""
function milky_way_orbit_parameters(V_potential)
    # Unit system: solar mass units where 1 M_sun = 1476.67 m and c = 1
    Msol2m = 1476.67           # 1 solar mass = 1476.67 m
    c_mks = 299792458.0        # Speed of light in m/s
    c = 1.0                    # Set c = 1 in our units
    
    # Convert standard units to solar mass units
    kpc_m = 3.0857e19          # 1 kpc in meters
    kpc = kpc_m / Msol2m       # 1 kpc in solar mass units
    km_s = 1000.0 / (Msol2m * c_mks)  # 1 km/s in solar mass units (with c=1)
    
    # Solar orbital parameters
    r_sun_kpc = 8.4            # Solar galactocentric distance in kpc
    r_sun = r_sun_kpc * kpc    # Convert to solar mass units
    
    # Calculate actual circular velocity from MW potential gradient
    x_sun = [r_sun, 0.0, 0.0]  # Solar position
    grad_mw = ∂(V_potential, x_sun)  # Gradient at solar position
    a_mw = norm(grad_mw)       # Acceleration magnitude
    v_circ = sqrt(a_mw * r_sun)  # Circular velocity: v = sqrt(a*r)
    
    # Calculate orbital period
    T_orbit = 2π * r_sun / v_circ
    
    println("=== Milky Way Orbital Parameters ===")
    println("Unit system: 1 M_sun = $(Msol2m) m, c = 1")
    println("Solar orbital radius: $(r_sun_kpc) kpc = $(r_sun) M_sun units")
    println("MW potential acceleration: $(a_mw) (c=1 units)")
    println("Calculated circular velocity: $(v_circ / km_s) km/s = $(v_circ) (c=1 units)")
    println("Target velocity (observational): 220.0 km/s")
    println("Orbital period: $(T_orbit) time units")
    
    return (r_sun, v_circ, T_orbit, kpc, km_s)
end

# Use the corrected Milky Way potential with origin at (0,0,0)
V = ΦMilkyWay(Double64, 0.0, 0.0, 0.0)  # Place galactic center at coordinate origin

# Get solar orbital parameters using the MW potential
r_sun, v_circ, T_orbit, kpc, km_s = milky_way_orbit_parameters(V)

# Solar mass only
M_sun = 1.0  # Solar mass in solar mass units
m = [M_sun]  # Just the Sun

# Construct the potential energy function
U = UConstructor(V, m, zeros(3))

# Initial conditions for circular orbit in x-y plane
# The MW potential includes origin offset, so we place Sun at simple coordinates
q0 = [[r_sun, 0.0, 0.0]]  # Position: Sun at distance r_sun from coordinate origin
p0 = [[0.0, m[1]*v_circ, 0.0]]  # Momentum: circular velocity in y-direction

system = pomin.Particles(m, q0, p0)

tspan = (0.0, T_orbit * 2.0)  # Simulate for 2 complete orbits to check stability
params = pomin.ParametersJulia(tspan)

sol = pomin.solve(system, params, U)

#   PLOTTING

# Extract trajectory data for plotting
n_particles = length(m)
d = 3  # 3D space
n_steps = length(sol.t)

# Extract positions for the Sun over time
x_traj = [sol.u[i][1] for i in 1:n_steps]  # Sun x-position
y_traj = [sol.u[i][2] for i in 1:n_steps]  # Sun y-position
z_traj = [sol.u[i][3] for i in 1:n_steps]  # Sun z-position

# Calculate orbital radius and velocity over time
orbital_radius = [sqrt(x_traj[i]^2 + y_traj[i]^2 + z_traj[i]^2) for i in 1:n_steps]
vx_traj = [sol.u[i][4]/m[1] for i in 1:n_steps]  # x-velocity
vy_traj = [sol.u[i][5]/m[1] for i in 1:n_steps]  # y-velocity
vz_traj = [sol.u[i][6]/m[1] for i in 1:n_steps]  # z-velocity
orbital_speed = [sqrt(vx_traj[i]^2 + vy_traj[i]^2 + vz_traj[i]^2) for i in 1:n_steps]

# Create orbital trajectory plot
println("Creating orbital trajectory plots...")

# Print trajectory statistics
println("\n=== Orbital Analysis ===")
println("Initial orbital radius: $(orbital_radius[1] / kpc) kpc")
println("Final orbital radius: $(orbital_radius[end] / kpc) kpc")
println("Radius variation: $((maximum(orbital_radius) - minimum(orbital_radius)) / kpc) kpc")
println("Initial orbital speed: $(orbital_speed[1] / km_s) km/s")
println("Final orbital speed: $(orbital_speed[end] / km_s) km/s")
println("Target speed: 220.0 km/s")

# 1. Orbital trajectory in x-y plane (convert to kpc for readability)
p1 = plot(x_traj ./ kpc, y_traj ./ kpc, 
          label="Solar Orbit", 
          linewidth=2, 
          color=:orange,
          title="Solar Orbit in Milky Way Potential",
          xlabel="x (kpc)", ylabel="y (kpc)",
          aspect_ratio=:equal)

# Mark initial position
scatter!(p1, [x_traj[1] / kpc], [y_traj[1] / kpc], 
         label="Start", 
         markersize=6, 
         color=:green, 
         markershape=:circle)

# Mark coordinate origin (not necessarily galactic center due to MW potential offset)
scatter!(p1, [0], [0], 
         label="Coordinate Origin", 
         markersize=8, 
         color=:black, 
         markershape=:star)

# 2. Orbital radius vs time
p2 = plot(sol.t ./ T_orbit, orbital_radius ./ kpc,
          label="Orbital Radius",
          linewidth=2,
          color=:blue,
          title="Orbital Radius vs Time",
          xlabel="Time (orbital periods)", ylabel="Radius (kpc)")

# Add expected radius line
hline!(p2, [r_sun / kpc], 
       label="Expected Radius (8.4 kpc)", 
       linestyle=:dash, 
       color=:red, 
       linewidth=2)

# 3. Orbital speed vs time
p3 = plot(sol.t ./ T_orbit, orbital_speed ./ km_s,
          label="Orbital Speed",
          linewidth=2,
          color=:green,
          title="Orbital Speed vs Time",
          xlabel="Time (orbital periods)", ylabel="Speed (km/s)")

# Add expected speed line
hline!(p3, [220.0], 
       label="Target Speed (220 km/s)", 
       linestyle=:dash, 
       color=:red, 
       linewidth=2)

# 4. 3D trajectory
p4 = plot(x_traj ./ kpc, y_traj ./ kpc, z_traj ./ kpc,
          label="Solar Orbit",
          linewidth=2,
          color=:purple,
          title="3D Solar Orbit",
          xlabel="x (kpc)", ylabel="y (kpc)", zlabel="z (kpc)")

# Combine plots
combined_plot = plot(p1, p2, p3, p4, layout=(2,2), size=(1000,800))

# Display the plot
display(combined_plot)

# Save the plot
savefig(combined_plot, "milky_way_potential_orbit.png")
println("Plot saved as 'milky_way_potential_orbit.png'")

# Final orbital statistics
println("\n=== Final Orbital Statistics ===")
println("Unit system: 1 M_sun = 1476.67 m, c = 1")
println("Simulation time: $(sol.t[end] / T_orbit) orbital periods")
println("Number of time steps: $n_steps")
println("Radius stability (max - min): $((maximum(orbital_radius) - minimum(orbital_radius)) / kpc) kpc")
println("Speed stability (max - min): $((maximum(orbital_speed) - minimum(orbital_speed)) / km_s) km/s")
println("Relative radius variation: $((maximum(orbital_radius) - minimum(orbital_radius)) / r_sun * 100)%")
println("Relative speed variation: $((maximum(orbital_speed) - minimum(orbital_speed)) / (220.0 * km_s) * 100)%")
println("Average orbital speed: $(sum(orbital_speed)/length(orbital_speed) / km_s) km/s")
println("Target speed: 220.0 km/s")

return sol, orbital_radius, orbital_speed
