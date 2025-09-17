include("../pomin.jl")
using .pomin
using LinearAlgebra, DoubleFloats, Plots

tpfl = Double64

#-----------------------------------------------------------------------
#   CONSTANTS
#-----------------------------------------------------------------------
cMKS = tpfl(299792458.0)
AUMKS = tpfl(149597870700.0)
lyMKS = tpfl(9460730472580800.0)
Msol2m = tpfl(1476.67)

c = one(tpfl)
AU = AUMKS / Msol2m
ly = lyMKS / Msol2m

#-----------------------------------------------------------------------
#   INITIAL DATA SETUP FOR SUN
#-----------------------------------------------------------------------

# Sun parameters (at origin)
msol = tpfl(1.0)
qsol = tpfl.([0.0, 0.0, 0.0])
psol = tpfl.([0.0, 0.0, 0.0])

Psol = pomin.setup_single_particle(msol, qsol, psol, tpfl)

#-----------------------------------------------------------------------
#   INITIAL DATA SETUP FOR JUPITER
#-----------------------------------------------------------------------

mjup = tpfl(0.000954)  # Jupiter mass in solar masses
qjup = tpfl.([3.99442362 * AU, 2.73345644 * AU, 1.07451704 * AU])  # Jupiter position in geometric units (from astropy, J2000)

vjup = tpfl.([-7.8875477321829800, 1.0175858402135900E+01, 4.5538873116399600])  # in km/s, from astropy, J2000

vjup *= 1000 / cMKS  # convert from km/s to units of c

γjup = one(tpfl) / sqrt(one(tpfl) - norm(vjup)^2)  # Lorentz factor
pjup = tpfl.(mjup * γjup * vjup)  # Jupiter momentum 

Pjup = pomin.setup_single_particle(mjup, qjup, pjup, tpfl)

#-----------------------------------------------------------------------
#   SETUP MAIN PARTICLE SYSTEM
#-----------------------------------------------------------------------

main_particle_system = pomin.merge_particle_systems(Psol)

#-----------------------------------------------------------------------
#   SETUP TEST PARTICLE SYSTEM
#-----------------------------------------------------------------------

test_particle_system = pomin.merge_particle_systems(Pjup)

#-----------------------------------------------------------------------
#   INTEGRATION PARAMETERS
#-----------------------------------------------------------------------
tspan = (tpfl(0), tpfl(1.6E+14)) # about 25 years (more than 2 orbits of Jupiter)
δ = tpfl(7.31E+8) # about one hour
tols = tpfl(1e-16)

params = pomin.ParametersJulia(tspan, integrator="Vern9",
    atol=tols, rtol=tols)

#-----------------------------------------------------------------------
#   RUN SOLVER
#-----------------------------------------------------------------------

println("Running solver...")
sol = pomin.solveT(main_particle_system, params, test_particle_system)

#   PLOTTING

# Extract trajectory data for plotting
n_particles = 2
d = 3  # 3D space
n_steps = length(sol.t)

# Extract positions for each particle over time
x1_traj = [sol.u[i][1] for i in 1:n_steps]  # Particle 1 x-position
y1_traj = [sol.u[i][2] for i in 1:n_steps]  # Particle 1 y-position
z1_traj = [sol.u[i][3] for i in 1:n_steps]  # Particle 1 z-position

x2_traj = [sol.u[i][7] for i in 1:n_steps]  # Particle 2 x-position
y2_traj = [sol.u[i][8] for i in 1:n_steps]  # Particle 2 y-position
z2_traj = [sol.u[i][9] for i in 1:n_steps]  # Particle 2 z-position

# Create orbital trajectory plot
println("Creating orbital trajectory plots...")

# 2D projection in x-y plane
p1 = plot(x1_traj, y1_traj,
    label="Sun",
    linewidth=2,
    color=:blue,
    title="Jupiter Orbit",
    xlabel="x", ylabel="y",
    aspect_ratio=:equal)
plot!(p1, x2_traj, y2_traj,
    label="Jupiter",
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
    label="Sun",
    markersize=8,
    color=:black,
    markershape=:star)

# 3D trajectory plot
p2 = plot(x1_traj, y1_traj, z1_traj,
    label="Sun",
    linewidth=2,
    color=:blue,
    title="3D Orbits",
    xlabel="x", ylabel="y", zlabel="z")
plot!(p2, x2_traj, y2_traj, z2_traj,
    label="Jupiter",
    linewidth=2,
    color=:red)

# Mark central mass at origin in 3D
scatter!(p2, [0], [0], [0],
    label="Sun",
    markersize=8,
    color=:black,
    markershape=:star)

# Time evolution of separation
separation = [sqrt((x1_traj[i] - x2_traj[i])^2 +
                   (y1_traj[i] - y2_traj[i])^2 +
                   (z1_traj[i] - z2_traj[i])^2) for i in 1:n_steps]

p3 = plot(sol.t, separation,
    label="Separation",
    linewidth=2,
    color=:green,
    title="Separation vs Time",
    xlabel="Time", ylabel="Separation")

# Combine plots
combined_plot = plot(p1, p2, p3, layout=(2, 2), size=(800, 600))

# Display the plot
display(combined_plot)

# Save the plot
savefig(combined_plot, "Jupiter_orbit.png")
println("Plot saved as 'Jupiter_orbit.png'")

# Print some orbital statistics
println("\n=== Orbital Statistics ===")
println("Integration time: $(sol.t[end]) time units")
println("Number of time steps: $n_steps")
println("Initial separation: $(separation[1])")
println("Final separation: $(separation[end])")
println("Max separation: $(maximum(separation))")
println("Min separation: $(minimum(separation))")
println("Average separation: $(sum(separation)/length(separation))")
