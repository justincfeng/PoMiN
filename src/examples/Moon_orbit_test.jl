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

# Moon coordinates generated using astropy with epoch J2000
qEarth = tpfl.([-1.8667826140E+07, 8.9633743993E+07, 3.8883291714E+07])         # geometric units
pEarth = tpfl.([-2.9839042080E-10, -5.0388889700E-11, -2.1846049145E-11])       # geometric units

PEarth = pomin.setup_single_particle(mEarth, qEarth, pEarth, tpfl)

#-----------------------------------------------------------------------
#   INITIAL DATA SETUP FOR MOON
#-----------------------------------------------------------------------

mMoon = 7.34767309E22 / 1.988416E30     # in units of solar masses

# Moon coordinates generated using astropy with epoch J2000
qMoon = tpfl.([-1.886529874442160E+07, 8.94531273942814E+07, 3.883175807801640E+07])    # geometric units
vMoon = tpfl.([-2.9141440121218000E+01, -5.6958177245812600, -2.4819741171379700]) ./ cMKS      # units of c
γMoon = one(tpfl) / sqrt(one(tpfl) - norm(vMoon)^2)  # Lorentz factor
pMoon = tpfl.(mMoon * γMoon * vMoon)  # Moon momentum 
   
PMoon = pomin.setup_single_particle(mMoon, qMoon, pMoon, tpfl)

#-----------------------------------------------------------------------
#   SETUP MAIN PARTICLE SYSTEM
#-----------------------------------------------------------------------

main_particle_system = pomin.merge_particle_systems(Psol, PEarth)

#-----------------------------------------------------------------------
#   SETUP TEST PARTICLE SYSTEM
#-----------------------------------------------------------------------

test_particle_system = pomin.merge_particle_systems(PMoon)

#-----------------------------------------------------------------------
#   INTEGRATION PARAMETERS
#-----------------------------------------------------------------------
tspan = (tpfl(0), tpfl(2.56E+13)) # about 4 years
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

n_steps = length(sol.t)

# Extract positions for each particle over time

# Earth
x1_traj = [sol.u[i][4] for i in 1:n_steps]  # Particle 1 x-position
y1_traj = [sol.u[i][5] for i in 1:n_steps]  # Particle 1 y-position
z1_traj = [sol.u[i][6] for i in 1:n_steps]  # Particle 1 z-position
# Moon
x2_traj = [sol.u[i][13] for i in 1:n_steps]  # Particle 2 x-position
y2_traj = [sol.u[i][14] for i in 1:n_steps]  # Particle 2 y-position
z2_traj = [sol.u[i][15] for i in 1:n_steps]  # Particle 2 z-position

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

# Display the plot
display(p3)

# Save the plot
savefig(p3, "Moon_separation.png")
println("Plot saved as 'Moon_separation.png'")

# Print some orbital statistics
println("\n=== Orbital Statistics ===")
println("Integration time: $(sol.t[end]) time units")
println("Number of time steps: $n_steps")
println("Initial separation: $(separation[1])")
println("Final separation: $(separation[end])")
println("Max separation: $(maximum(separation))")
println("Min separation: $(minimum(separation))")
println("Average separation: $(sum(separation)/length(separation))")
