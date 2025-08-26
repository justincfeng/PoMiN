#-----------------------------------------------------------------------
#
#   SUN-JUPITER ORBIT EXAMPLE WITH MANUAL INITIAL DATA SETUP
#
#-----------------------------------------------------------------------

#!/usr/bin/env julia

include("../pomin.jl")
using .pomin
using LinearAlgebra, Plots

#-----------------------------------------------------------------------
#   PARAMETERS
#-----------------------------------------------------------------------

m1 = 1.0        # Mass of Sun (normalized)
m2 = 0.001      # Mass of Jupiter (approx. 1/1000 of Sun's mass)
r = 5.2e9       # Orbital radius (5.2 AU in meters, scaled for simulation)
G = 1.0         # Gravitational constant (normalized for simulation)

# Calculate orbital period using Kepler's third law: T = 2π√(a³/GM)
M = m1 + m2                    # Total mass
a = r                          # For circular orbit, semi-major axis = radius
T = 2π * sqrt(a^3 / (G * M))   # Orbital period
n_orbits = 1.0                 # Number of orbits to simulate

# Time span
t_start = 0.0
t_end   = n_orbits * T

# Set up integration parameters
params = ParametersJulia((t_start, t_end))

#-----------------------------------------------------------------------
#   SET UP INITIAL DATA MANUALLY
#-----------------------------------------------------------------------

# Setup Sun as the first particle at the origin with zero momentum
sun_pos = [0.0, 0.0, 0.0]  # Position at origin
sun_mom = [0.0, 0.0, 0.0]  # Zero initial momentum
particles = pomin.setup_single_particle(m1, sun_pos, sun_mom)

# Calculate Jupiter's initial position and velocity for circular orbit
# For simplicity, place Jupiter at (r * m1/(m1+m2), 0, 0) adjusted for center of mass
jupiter_x = r * m1 / (m1 + m2)
jupiter_pos = [jupiter_x, 0.0, 0.0]
# Calculate orbital velocity for circular orbit: v = sqrt(G * m1^2 / (r * (m1 + m2))) for Jupiter relative to center of mass
# Using the reduced mass concept for two-body circular orbit
v_jupiter = sqrt(G * m1^2 / (r * (m1 + m2)))  # Velocity of Jupiter relative to center of mass
jupiter_mom = [0.0, m2 * v_jupiter, 0.0]  # Momentum = mass * velocity, in y-direction for circular orbit

# Add Jupiter to the system
particles = pomin.add_particle(particles, m2, jupiter_pos, jupiter_mom)

# Add Earth as a test particle at 1 AU from Sun
m_earth = 3.0e-6  # Small mass compared to Sun, acting as a test particle for visibility
r_earth = 1.496e9  # 1 AU in meters, scaled similarly to Jupiter's radius
# Position Earth at 1 AU, adjusted for center of mass (approximately at Sun's position since Sun is much heavier)
earth_x = r_earth * m1 / (m1 + m2)  # Very close to r_earth since m2 << m1
# Place Earth at 90 degrees ahead of Jupiter for visibility (along y-axis for simplicity)
earth_pos = [0.0, earth_x, 0.0]
# Calculate orbital velocity for Earth: v = sqrt(G * m1 / r_earth) since m_earth is negligible
v_earth = sqrt(G * m1 / r_earth)
earth_mom = [m_earth * v_earth, 0.0, 0.0]  # Momentum in x-direction for circular orbit

# Add Earth to the system
testparticles = Particles([m_earth], [[earth_x * cos(60 * pi/180), earth_x * sin(60 * pi/180), 0.0]], [[-v_earth * sin(60 * pi/180), v_earth * cos(60 * pi/180), 0.0]])

# Extract momentum and positions
total_momentum = particles.p[1] + particles.p[2]  # Only Sun and Jupiter for main system
com_position = (particles.m[1] * particles.q[1] + particles.m[2] * particles.q[2]) / (particles.m[1] + particles.m[2])

# Adjust Sun's momentum to ensure total momentum is zero (center of mass frame)
# Total momentum should be zero, so Sun's momentum = -(Jupiter's momentum)
particles.p[1] = -particles.p[2]

#-----------------------------------------------------------------------
#   SOLVE PROBLEM
#-----------------------------------------------------------------------

@time solution = pomin.solveT(particles, params, testparticles)

n_steps        = length(solution.t)
times          = solution.t

x1  = [solution.u[i][1] for i in 1:n_steps]  # x-position of Sun
y1  = [solution.u[i][2] for i in 1:n_steps]  # y-position of Sun
x2  = [solution.u[i][4] for i in 1:n_steps]  # x-position of Jupiter  
y2  = [solution.u[i][5] for i in 1:n_steps]  # y-position of Jupiter
x3  = [solution.u[i][13] for i in 1:n_steps]  # x-position of Earth  
y3  = [solution.u[i][14] for i in 1:n_steps]  # y-position of Earth

# Extract momenta for energy calculation
px1 = [solution.u[i][7] for i in 1:n_steps]  # x-momentum of Sun
py1 = [solution.u[i][8] for i in 1:n_steps]  # y-momentum of Sun
px2 = [solution.u[i][10] for i in 1:n_steps] # x-momentum of Jupiter
py2 = [solution.u[i][11] for i in 1:n_steps] # y-momentum of Jupiter
px3 = [solution.u[i][16] for i in 1:n_steps] # x-momentum of Earth
py3 = [solution.u[i][17] for i in 1:n_steps] # y-momentum of Earth

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
ET0 = HamPM.HT(vcat(particles.q..., particles.p...), particles.m, vcat(testparticles.q..., testparticles.p...), testparticles.m, 3)
println("Initial energy for test system: ", ET0)

initial_energy = energies[1]
final_energy = energies[end]
energy_drift = abs(final_energy - initial_energy) / abs(initial_energy)

# Plot 1: All orbital trajectories on same xy plot
p1 = plot(aspect_ratio=:equal,
    title="Sun-Jupiter-Earth System Orbit (Manual Setup)",
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

savefig(p1, "pomin_sun_jupiter_earth_manual_orbit.pdf")

# Plot 2: Energy conservation
p2 = plot(times, energies, 
    title="Energy Conservation (Sun-Jupiter-Earth System)",
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
savefig(p2, "pomin_sun_jupiter_earth_manual_energy_conservation.pdf")

# Plot 3: Separation between Sun-Jupiter and Sun-Earth
separation12 = sqrt.((x2 .- x1).^2 + (y2 .- y1).^2)
separation13 = sqrt.((x3 .- x1).^2 + (y3 .- y1).^2)
p3 = plot(times, separation12, 
    title="Separation in Sun-Jupiter-Earth System",
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
savefig(p3, "pomin_sun_jupiter_earth_manual_separation.pdf")

nothing
