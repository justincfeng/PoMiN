

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
    solar_circular_orbit_conditions(potential_func)
    
Returns initial conditions for the Sun's circular orbit in the Milky Way.
Coordinates are in solar mass units to match the potential function.
"""
function solar_circular_orbit_conditions(potential_func)
    # Unit conversions (geometric units where G=c=1)
    kpc = 5.224206385e15  # kpc in units of solar masses
    km_s = 3.335640952e-6  # km/s in units of c
    
    # Solar position relative to galactic center (in solar mass units)
    # The potential function adds origin offset, so we need to account for that
    origin_x = -1.708859462494220E+17  # Sun's position in galactic frame
    origin_y = 0.0
    origin_z = 4.346342845091530E+14
    
    # Solar galactocentric distance in kpc, then convert to solar mass units
    r_sun_kpc = 8.4  # kpc
    r_sun_sim = r_sun_kpc * kpc  # Convert to solar mass units
    
    # Use observational circular velocity (known to work)
    v_circ = 220.0 * km_s  # Observational circular velocity
    
    # Calculate orbital period
    T_orbit = 2π * r_sun_sim / v_circ
    
    return (r_sun_sim, v_circ, T_orbit)
end #-------------------------------------------------------------------

# Test with simple central mass potential first
V_test = ΦCentralMass(zeros(3), Double64, 2.325e7*3283.0, 1.0)  # Total galactic mass at center
# V = ΦMilkyWay(Double64)  # Comment out for now

# Get solar orbital parameters using the test potential
r_sun, v_circ, T_orbit = solar_circular_orbit_conditions(V_test)

# Solar mass only
M_sun = 1.0  # Solar mass in solar mass units
m = [M_sun]  # Just the Sun

# Construct the potential energy function
U = UConstructor(V_test, m, zeros(3))

# Place Sun at 8.4 kpc from galactic center (which is at origin due to potential offset)
# Simple circular orbit in x-z plane
q0 = [[r_sun, 0.0, 0.0]]  # Sun at (8.4 kpc, 0, 0) from galactic center

# For circular orbit: velocity perpendicular to radius
# Position [r, 0, 0] -> velocity [0, 0, v] for counterclockwise orbit in x-z plane
p0 = [[0.0, 0.0, m[1]*v_circ]]  # Circular velocity in +z direction

system = pomin.Particles(m, q0, p0)

tspan = (0.0, T_orbit * 1.1)  # Simulate for slightly more than one complete orbit
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

# Create orbital trajectory plot
println("Creating orbital trajectory plots...")

# Print trajectory statistics to understand the motion
println("Trajectory analysis:")
println("X range: $(minimum(x_traj)) to $(maximum(x_traj))")
println("Y range: $(minimum(y_traj)) to $(maximum(y_traj))")
println("Z range: $(minimum(z_traj)) to $(maximum(z_traj))")

# Center the plot on the Sun's orbit (subtract mean position)
x_mean = sum(x_traj) / length(x_traj)
y_mean = sum(y_traj) / length(y_traj)
z_mean = sum(z_traj) / length(z_traj)

x_centered = x_traj .- x_mean
y_centered = y_traj .- y_mean
z_centered = z_traj .- z_mean

println("Centered ranges:")
println("X centered: $(minimum(x_centered)) to $(maximum(x_centered))")
println("Y centered: $(minimum(y_centered)) to $(maximum(y_centered))")
println("Z centered: $(minimum(z_centered)) to $(maximum(z_centered))")

# 2D projection in galactic (x-z) plane using centered coordinates
p1 = plot(x_centered, z_centered, 
          label="Sun", 
          linewidth=2, 
          color=:orange,
          title="Solar Orbit in Galactic Plane (centered)",
          xlabel="x - x_mean (solar masses)", ylabel="z - z_mean (solar masses)",
          aspect_ratio=:equal)

# Mark initial position
scatter!(p1, [x_centered[1]], [z_centered[1]], 
         label="Sun Start", 
         markersize=6, 
         color=:orange, 
         markershape=:circle)

# Mark center of orbit
scatter!(p1, [0], [0], 
         label="Orbit Center", 
         markersize=8, 
         color=:black, 
         markershape=:star)

# 3D trajectory plot using centered coordinates
p2 = plot(x_centered, y_centered, z_centered,
          label="Sun",
          linewidth=2,
          color=:orange,
          title="3D Solar Orbit (centered)",
          xlabel="x - x_mean", ylabel="y - y_mean", zlabel="z - z_mean")

# Time evolution of orbital radius (distance from mean position)
orbital_radius = [sqrt(x_centered[i]^2 + y_centered[i]^2 + z_centered[i]^2) for i in 1:n_steps]

p3 = plot(sol.t, orbital_radius,
          label="Orbital Radius",
          linewidth=2,
          color=:green,
          title="Orbital Radius vs Time",
          xlabel="Time (geometric units)", ylabel="Radius from orbit center")

# Additional plot: x-y plane motion (centered)
p4 = plot(x_centered, y_centered,
          label="Sun",
          linewidth=2,
          color=:red,
          title="Motion in x-y plane (centered)",
          xlabel="x - x_mean", ylabel="y - y_mean",
          aspect_ratio=:equal)

# Mark orbit center
scatter!(p4, [0], [0],
         label="Orbit Center",
         markersize=8,
         color=:black,
         markershape=:star)

# Combine plots
combined_plot = plot(p1, p2, p3, p4, layout=(2,2), size=(1000,800))

# Display the plot
display(combined_plot)

# Save the plot
savefig(combined_plot, "solar_orbit_milky_way_potential.png")
println("Plot saved as 'solar_orbit_milky_way_potential.png'")

# Print some orbital statistics
println("\n=== Solar Orbital Statistics ===")
println("Integration time: $(sol.t[end]) time units")
println("Number of time steps: $n_steps")
println("Initial orbital radius: $(orbital_radius[1])")
println("Final orbital radius: $(orbital_radius[end])")
println("Max orbital radius: $(maximum(orbital_radius))")
println("Min orbital radius: $(minimum(orbital_radius))")
println("Average orbital radius: $(sum(orbital_radius)/length(orbital_radius))")
println("Circular velocity: $v_circ")
println("Orbital period: $T_orbit")
