#-----------------------------------------------------------------------
#
#   SUN-JUPITER ORBIT EXAMPLE WITH MANUAL INITIAL DATA SETUP
#
#-----------------------------------------------------------------------

#!/usr/bin/env julia

include("../pomin.jl")
using .pomin
using LinearAlgebra, Plots, Printf

#-----------------------------------------------------------------------
#   PARAMETERS
#-----------------------------------------------------------------------

m1 = 1.0        # Mass of Sun (normalized)
m2 = 0.001      # Mass of Jupiter (approx. 1/1000 of Sun's mass)
m3 = 1.0       # Mass of Earth (small test particle mass)
r = 5.2e8         # Orbital radius (5.2 AU in normalized units)
G = 1.0         # Gravitational constant (normalized for simulation)

# Calculate orbital period using Kepler's third law: T = 2π√(a³/GM)
M = m1 + m2                    # Total mass
a = r                          # For circular orbit, semi-major axis = radius
T = 2π * sqrt(a^3 / (G * M))   # Orbital period
n_orbits = 2.0                 # Number of orbits to simulate

println("Orbital parameters:")
println("  Jupiter orbital radius: $r")
println("  Total mass: $M")
println("  Orbital period: $T")
println("  Simulation time: $(n_orbits * T)")

# Time span
t_start = 0.0
t_end   = n_orbits * T

# Set up integration parameters
params = ParametersJulia((t_start, t_end))

#-----------------------------------------------------------------------
#   SET UP INITIAL DATA MANUALLY
#-----------------------------------------------------------------------

# Calculate positions and velocities for circular orbit around center of mass
# Jupiter position relative to center of mass
jupiter_x = r * m1 / (m1 + m2)
# Sun position relative to center of mass (opposite direction)
sun_x = -r * m2 / (m1 + m2)

# Circular orbital velocity for reduced mass system
v_orbit = sqrt(G * (m1 + m2) / r)

# Sun
sun_pos = [sun_x, 0.0, 0.0]  # Position offset from origin
sun_mom = [0.0, m1 * v_orbit * m2 / (m1 + m2), 0.0]  # Counter-momentum to Jupiter

# Jupiter  
jupiter_pos = [jupiter_x, 0.0, 0.0]
jupiter_mom = [0.0, -m2 * v_orbit * m1 / (m1 + m2), 0.0]  # Opposite momentum to Sun

# Initialize system with Sun
particles = pomin.setup_single_particle(m1, sun_pos, sun_mom)

# Add Jupiter to the system
particles = pomin.add_particle(particles, m2, jupiter_pos, jupiter_mom)

# Earth (test particle) - using smaller radius for better visibility
m_earth = m3  # Use the defined mass for Earth
r_earth = 1e8  # Smaller radius for more visible motion
earth_pos = [r_earth, 0.0, 0.0]  # Place Earth at smaller radius
v_earth = sqrt(G * m1 / r_earth)  # Circular orbital velocity
earth_mom = [0.0, m_earth * v_earth, 0.0]  # Momentum in y-direction for circular orbit

println("Earth orbital setup:")
println("  Earth orbital radius: $r_earth")
println("  Earth orbital velocity: $v_earth")
println("  Earth momentum magnitude: $(m_earth * v_earth)")

# Add Earth to the system
testparticles = Particles([m_earth], [earth_pos], [earth_mom])

# Check momentum and center of mass
total_momentum = particles.p[1] + particles.p[2]  # Should be zero for center of mass frame
com_position = (particles.m[1] * particles.q[1] + particles.m[2] * particles.q[2]) / (particles.m[1] + particles.m[2])

println("System check:")
println("  Total momentum: $total_momentum")
println("  Center of mass position: $com_position")

#-----------------------------------------------------------------------
#   SOLVE PROBLEM
#-----------------------------------------------------------------------

println("\nInitial conditions:")
println("  Sun position: $(particles.q[1])")
println("  Sun momentum: $(particles.p[1])")
println("  Jupiter position: $(particles.q[2])")  
println("  Jupiter momentum: $(particles.p[2])")
println("  Earth position: $(testparticles.q[1])")
println("  Earth momentum: $(testparticles.p[1])")



@time solution = pomin.solve(particles, params; testparticles=testparticles)

n_steps        = length(solution.t)
times          = solution.t

# Phase space layout: [q1, q2, qT, p1, p2, pT] where each q,p is 3D
# Sun (particle 1): positions [1:3], momenta [10:12]
# Jupiter (particle 2): positions [4:6], momenta [13:15] 
# Earth (test particle): positions [7:9], momenta [16:18]

x1  = [solution.u[i][1] for i in 1:n_steps]   # x-position of Sun
y1  = [solution.u[i][2] for i in 1:n_steps]   # y-position of Sun
x2  = [solution.u[i][4] for i in 1:n_steps]   # x-position of Jupiter  
y2  = [solution.u[i][5] for i in 1:n_steps]   # y-position of Jupiter
x3  = [solution.u[i][13] for i in 1:n_steps]   # x-position of Earth  
y3  = [solution.u[i][14] for i in 1:n_steps]   # y-position of Earth

# Extract momenta for energy calculation
px1 = [solution.u[i][6] for i in 1:n_steps]  # x-momentum of Sun
py1 = [solution.u[i][7] for i in 1:n_steps]  # y-momentum of Sun
px2 = [solution.u[i][9] for i in 1:n_steps]  # x-momentum of Jupiter
py2 = [solution.u[i][10] for i in 1:n_steps]  # y-momentum of Jupiter
px3 = [solution.u[i][16] for i in 1:n_steps]  # x-momentum of Earth
py3 = [solution.u[i][17] for i in 1:n_steps]  # y-momentum of Earth

# Calculate energy at each time step (only for main system: Sun and Jupiter)
energies = Float64[]
for i in 1:n_steps
    z = vcat(particles.q[1], particles.q[2], particles.p[1], particles.p[2])  # Extract only Sun and Jupiter data (first two bodies, 3D position + momentum per body)
    masses = particles.m  # Only main system masses (Sun and Jupiter)
    H = HamPM.H(vcat(particles.q..., particles.p...), masses, 3)  # Dimension 3 for 2D problem with z-component=0
    push!(energies, H)
end

println("Main system state vector length: ", length(vcat(particles.q..., particles.p...)))
println("Test system state vector length: ", length(vcat(testparticles.q..., testparticles.p...)))
E0 = HamPM.H(vcat(particles.q..., particles.p...), particles.m, 3)
println("Initial energy for main system: ", E0)
ET0 = HamPM.HT(vcat(testparticles.q..., testparticles.p...), testparticles.m, vcat(particles.q..., particles.p...), particles.m, 3)
println("Initial energy for test system: ", ET0)

initial_energy = energies[1]
final_energy = energies[end]
energy_drift = abs(final_energy - initial_energy) / abs(initial_energy)

# Check orbital motion by examining positions at different times
println("\nOrbital motion check:")
n_check = min(5, n_steps)
for i in 1:n_check
    idx = Int(round((i-1) * (n_steps-1) / (n_check-1))) + 1
    t = times[idx]
    println("  t=$(round(t, digits=2)): Jupiter=($(@sprintf("%.3f", x2[idx])), $(@sprintf("%.3f", y2[idx]))), Earth=($(@sprintf("%.6f", x3[idx])), $(@sprintf("%.6f", y3[idx])))")
end

# Check if Earth is actually moving (with higher precision)
earth_x_change = abs(x3[end] - x3[1])
earth_y_change = abs(y3[end] - y3[1])
println("\nEarth motion analysis:")
println("  Initial Earth position: ($(@sprintf("%.10f", x3[1])), $(@sprintf("%.10f", y3[1])))")
println("  Final Earth position: ($(@sprintf("%.10f", x3[end])), $(@sprintf("%.10f", y3[end])))")
println("  Earth x-displacement: $(@sprintf("%.2e", earth_x_change))")
println("  Earth y-displacement: $(@sprintf("%.2e", earth_y_change))")
println("  Earth is moving: $(earth_x_change > 1e-6 || earth_y_change > 1e-6)")

# Calculate Jupiter's distance from Sun over time
jupiter_distances = [sqrt(x2[i]^2 + y2[i]^2) for i in 1:n_steps]
earth_distances = [sqrt(x3[i]^2 + y3[i]^2) for i in 1:n_steps]

println("\nOrbital radii:")
println("  Jupiter: min=$(round(minimum(jupiter_distances), digits=3)), max=$(round(maximum(jupiter_distances), digits=3)), avg=$(round(sum(jupiter_distances)/length(jupiter_distances), digits=3))")
println("  Earth: min=$(round(minimum(earth_distances), digits=3)), max=$(round(maximum(earth_distances), digits=3)), avg=$(round(sum(earth_distances)/length(earth_distances), digits=3))")

# Plot 1: All orbital trajectories on same xy plot
p1 = plot(aspect_ratio=:equal,
    title="Sun-Jupiter-(test)Earth System Orbit (Manual Setup)",
    xlabel="x position", 
    ylabel="y position",
    legend=:topright,
    grid=true,
    gridwidth=1,
    gridcolor=:gray,
    gridalpha=0.5,
    background_color_inside=:white,
    foreground_color_grid=:gray,
    size=(800,600)
)

# Plot Sun's trajectory
plot!(p1, x1, y1, label="Sun", color=:goldenrod, linewidth=1.5)
# Plot Jupiter's trajectory
plot!(p1, x2, y2, label="Jupiter", color=:orangered, linewidth=1.5)
# Plot Earth's trajectory
plot!(p1, x3, y3, label="Earth", color=:blue, linewidth=1.5)

# Mark initial positions
scatter!(p1, [x1[1]], [y1[1]], label=false, color=:goldenrod, markersize=5, marker=:circle)
scatter!(p1, [x2[1]], [y2[1]], label=false, color=:orangered, markersize=5, marker=:circle)
scatter!(p1, [x3[1]], [y3[1]], label=false, color=:blue, markersize=5, marker=:circle)

savefig(p1, "pomin_sun_jupiter_test-earth_orbit.pdf")

# Plot 2: Energy conservation
p2 = plot(times, energies, 
    title="Energy Conservation (Sun-Jupiter-(test)_Earth System)",
    xlabel="Time",
    ylabel="Hamiltonian Energy",
    label="Energy",
    color=:red,
    linewidth=1.5,
    grid=true,
    gridwidth=1,
    gridcolor=:gray,
    gridalpha=0.5,
    background_color_inside=:white,
    foreground_color_grid=:gray,
    size=(800,400)
)
savefig(p2, "pomin_sun_jupiter_test-earth_energy_conservation.pdf")

# Plot 3: Separation between Sun-Jupiter and Sun-Earth
separation12 = sqrt.((x2 .- x1).^2 + (y2 .- y1).^2)
separation13 = sqrt.((x3 .- x1).^2 + (y3 .- y1).^2)
p3 = plot(times, separation12, 
    title="Separation in Sun-Jupiter-(test)_Earth System",
    xlabel="Time",
    ylabel="Separation",
    label="Sun-Jupiter",
    color=:orangered,
    linewidth=1.5,
    grid=true,
    gridwidth=1,
    gridcolor=:gray,
    gridalpha=0.5,
    background_color_inside=:white,
    foreground_color_grid=:gray,
    size=(800,400)
)
plot!(p3, times, separation13, label="Sun-Earth", color=:blue, linewidth=1.5)
savefig(p3, "pomin_sun_jupiter_test-earth_separation.pdf")

nothing
