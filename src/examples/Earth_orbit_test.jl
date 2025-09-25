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
#   INITIAL DATA SETUP FOR EARTH
#-----------------------------------------------------------------------

mEarth = tpfl(3.0033693739E-06)

# Earth coordinates generated using astropy with epoch J2000
qEarth = tpfl.([-1.866782584728260E+07, 8.963374397323200E+07, 3.888329146643650E+07])      # geometric units
vEarth = tpfl.([-9.935188906586840E-05, -1.677745335408000E-05, -7.273846929131220E-06])    # units of c
γEarth = one(tpfl) / sqrt(one(tpfl) - norm(vEarth)^2)
pEarth = mEarth * γEarth * vEarth

PEarth = pomin.setup_single_particle(mEarth, qEarth, pEarth, tpfl)

#-----------------------------------------------------------------------
#   SETUP MAIN PARTICLE SYSTEM
#-----------------------------------------------------------------------

main_particle_system = pomin.merge_particle_systems(Psol)

#-----------------------------------------------------------------------
#   SETUP TEST PARTICLE SYSTEM
#-----------------------------------------------------------------------

test_particle_system = pomin.merge_particle_systems(PEarth)

#-----------------------------------------------------------------------
#   INTEGRATION PARAMETERS
#-----------------------------------------------------------------------
tspan = (tpfl(0), tpfl(2.56E+14)) # about 40 years
δ = tpfl(7.31E+8) # about one hour
tols = tpfl(1e-16)

params = pomin.ParametersJulia(tspan, integrator="Vern9",
    atol=tols, rtol=tols)

#-----------------------------------------------------------------------
#   RUN SOLVER
#-----------------------------------------------------------------------

println("Running solver...")
sol = pomin.solve(main_particle_system, params; testparticles=test_particle_system)

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
    title="Earth Orbit",
    xlabel="x", ylabel="y",
    aspect_ratio=:equal)
plot!(p1, x2_traj, y2_traj,
    label="Earth",
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
    label="Earth",
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
savefig(combined_plot, "Earth_orbit.png")
println("Plot saved as 'Earth_orbit.png'")

# Print some orbital statistics
println("\n=== Orbital Statistics ===")
println("Integration time: $(sol.t[end]) time units")
println("Number of time steps: $n_steps")
println("Initial separation: $(separation[1])")
println("Final separation: $(separation[end])")
println("Max separation: $(maximum(separation))")
println("Min separation: $(minimum(separation))")
println("Average separation: $(sum(separation)/length(separation))")
