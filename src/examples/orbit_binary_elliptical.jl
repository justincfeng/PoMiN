#-----------------------------------------------------------------------
#
#   ELLIPTICAL ORBIT EXAMPLE
#
#-----------------------------------------------------------------------

#!/usr/bin/env julia

include("../pomin.jl")
using .pomin
using LinearAlgebra, Plots

#-----------------------------------------------------------------------
#   PARAMETERS
#-----------------------------------------------------------------------

m1 = 1.0        # Mass of first particle
m2 = 1.0        # Mass of second particle  
a = 1e6         # Semi-major axis
e = 0.7         # Eccentricity (0 < e < 1 for elliptical orbits)
G = 1.0         # Gravitational constant

# Calculate orbital period using Kepler's third law: T = 2π√(a³/GM)
M = m1 + m2                    # Total mass
T = 2π * sqrt(a^3 / (G * M))   # Orbital period
n_orbits = 0.5                 # Number of orbits to simulate

# Time span
t_start = 0.0
t_end   = n_orbits * T

# Set up integration parameters
params = ParametersJulia((t_start, t_end))

#-----------------------------------------------------------------------
#   SET UP INITIAL DATA
#-----------------------------------------------------------------------
particles = setup_elliptical_orbit(m1, m2, a, e; G=G)

# Extract momentum and positions
total_momentum = particles.p[1] + particles.p[2]
com_position = (particles.m[1] * particles.q[1] + particles.m[2] * particles.q[2]) / (particles.m[1] + particles.m[2])

#-----------------------------------------------------------------------
#   SOLVE PROBLEM
#-----------------------------------------------------------------------

@time solution = solve(particles, params)

n_steps        = length(solution.t)
times          = solution.t

x1  = [solution.u[i][1] for i in 1:n_steps]  # x-position of particle 1
y1  = [solution.u[i][2] for i in 1:n_steps]  # y-position of particle 1
x2  = [solution.u[i][4] for i in 1:n_steps]  # x-position of particle 2  
y2  = [solution.u[i][5] for i in 1:n_steps]  # y-position of particle 2

# Extract momenta for energy calculation
px1 = [solution.u[i][7] for i in 1:n_steps]  # x-momentum of particle 1
py1 = [solution.u[i][8] for i in 1:n_steps]  # y-momentum of particle 1
px2 = [solution.u[i][10] for i in 1:n_steps] # x-momentum of particle 2
py2 = [solution.u[i][11] for i in 1:n_steps] # y-momentum of particle 2

# Calculate energy at each time step
energies = Float64[]
for i in 1:n_steps
    z = solution.u[i]
    masses = particles.m
    H = HamPM.H(z, masses, 3)
    push!(energies, H)
end

initial_energy = energies[1]
final_energy = energies[end]
energy_drift = abs(final_energy - initial_energy) / abs(initial_energy)

# Plot 1: Both orbital trajectories on same xy plot
p1 = plot(aspect_ratio=:equal,
    title="Binary System Elliptical Orbits (e=$e)",
    xlabel="x position", 
    ylabel="y position",
    legend=:topright,
    grid=true,
    gridwidth=1,
    gridcolor=:lightgray)

# Plot both particle trajectories
plot!(p1, x1, y1, 
    label="Particle 1 (m=$m1)", 
    linewidth=3, 
    color=:blue,
    alpha=0.8)

plot!(p1, x2, y2, 
    label="Particle 2 (m=$m2)", 
    linewidth=3, 
    color=:red,
    alpha=0.8)

# Mark initial positions with larger markers (at apoapsis)
scatter!(p1, [x1[1]], [y1[1]], 
    color=:blue, 
    markersize=4, 
    markerstroke=2,
    markerstrokecolor=:darkblue,
    label="Start 1 (apoapsis)")
scatter!(p1, [x2[1]], [y2[1]], 
    color=:red, 
    markersize=4,
    markerstroke=2, 
    markerstrokecolor=:darkred,
    label="Start 2 (apoapsis)")

# Mark final positions
scatter!(p1, [x1[end]], [y1[end]], 
    color=:lightblue, 
    markersize=4,
    marker=:square,
    label="End 1")
scatter!(p1, [x2[end]], [y2[end]], 
    color=:pink, 
    markersize=4,
    marker=:square,
    label="End 2")

# Mark center of mass with prominent marker
scatter!(p1, [0], [0], 
    color=:black, 
    marker=:x, 
    markersize=6,
    markerstrokewidth=2,
    label="Center of Mass")

# Add ellipse showing the expected orbital shape for particle 1
θ = 0:0.01:2π
# Semi-minor axis
b = a * sqrt(1 - e^2)
# Distance from COM to particle 1 focus
c = a * e * m2 / (m1 + m2)
# Ellipse centered at particle 1's focus
ellipse_x = (a * m2 / (m1 + m2)) * cos.(θ) .- c
ellipse_y = (b * m2 / (m1 + m2)) * sin.(θ)
plot!(p1, ellipse_x, ellipse_y,
    color=:gray,
    linestyle=:dash,
    linewidth=1,
    alpha=0.5,
    label="Expected ellipse")

# Plot 2: Energy conservation
p2 = plot(times, energies, 
    label="Total Energy", 
    linewidth=2, 
    color=:green,
    title="Energy Conservation",
    xlabel="Time", 
    ylabel="Energy",
    aspect_ratio=:auto)

# Plot 3: Separation distance
separations = [sqrt((x1[i]-x2[i])^2 + (y1[i]-y2[i])^2) for i in 1:n_steps]
p3 = plot(times, separations, 
    label="Separation", 
    linewidth=2, 
    color=:purple,
    title="Particle Separation",
    xlabel="Time", 
    ylabel="Distance",
    aspect_ratio=:auto)

# Add expected periapsis and apoapsis lines
r_periapsis = a * (1 - e)
r_apoapsis = a * (1 + e)
hline!(p3, [r_periapsis], 
    label="Periapsis ($r_periapsis)", 
    linestyle=:dash, 
    color=:red)
hline!(p3, [r_apoapsis], 
    label="Apoapsis ($r_apoapsis)", 
    linestyle=:dash, 
    color=:blue)

# Display and save
savefig(p1, "pomin_elliptical_orbit.pdf")
savefig(p2, "pomin_elliptical_energy_conservation.pdf")
savefig(p3, "pomin_elliptical_separation.pdf")

# Print orbital parameters
println("Orbital Parameters:")
println("  Semi-major axis (a): $a")
println("  Eccentricity (e): $e")
println("  Periapsis distance: $r_periapsis")
println("  Apoapsis distance: $r_apoapsis")
println("  Orbital period (T): $T")
println("  Simulation time: $(n_orbits) orbits = $(t_end)")
println("  Energy drift: $(energy_drift*100)%")

nothing
