#!/usr/bin/env julia

#-----------------------------------------------------------------------
#
#   TEST SCRIPT FOR solveT WITH CENTRAL POTENTIAL
#   Two equal mass particles on opposite sides
#   One as main particle, one as test particle
#
#-----------------------------------------------------------------------

include("../pomin.jl")
using .pomin
using LinearAlgebra, Printf, Plots, Statistics
using DoubleFloats

tpfl = Double64

# Include external potential functionality
include("../core/physics/external_potentials/external.jl")

#-----------------------------------------------------------------------
#   CONSTANTS AND PARAMETERS
#-----------------------------------------------------------------------

# Massive interacting particle
M_massive = tpfl(1.0)   # Massive particle mass

# Test particle mass
m_test = tpfl(1e-6)     # Small test mass

# Initial positions
Rs = 2 * M_massive      # Schwarzschild radius
r0 = tpfl(1e6) * Rs     # ~1 million Schwarzschild radii
q_main = [tpfl(0.0), tpfl(0.0), tpfl(0.0)]  # Massive particle at origin
q_test = [r0, tpfl(0.0), tpfl(0.0)]         # Test particle at distance r0

# Circular orbital velocity for test particle
v_circular = sqrt(M_massive / r0)  # Circular velocity around massive particle
p_main = [tpfl(0.0), tpfl(0.0), tpfl(0.0)]                    # Massive particle at rest
p_test = [tpfl(0.0), m_test * v_circular, tpfl(0.0)]          # Test particle with circular velocity

# Massive interacting particle system
PintMain = pomin.setup_single_particle(M_massive, q_main, p_main, tpfl)

# Test particle system (feels forces but doesn't exert them)
PtestTest = pomin.setup_single_particle(m_test, q_test, p_test, tpfl)

# One orbital period for comparison
T_orbit = 2*π*sqrt(r0^3/M_massive)
tspan = (tpfl(0.0), T_orbit)
tols = tpfl(1e-12)

params = pomin.ParametersJulia(tspan, integrator="Vern9", 
                               atol=tols, rtol=tols)

sol_central = pomin.solveT(PintMain, params, PtestTest)

# Final positions
Z_central = sol_central(T_orbit)
q_main_final = Z_central[1:3]
q_test_final = Z_central[7:9]

#-----------------------------------------------------------------------
#   PLOTTING
#-----------------------------------------------------------------------

# Extract trajectory data
t_points = range(0, T_orbit, length=1000)
traj_central = [sol_central(t) for t in t_points]

# Main particle trajectories
x_main = [traj[1] for traj in traj_central]
y_main = [traj[2] for traj in traj_central]

# Test particle trajectories
x_test = [traj[7] for traj in traj_central]
y_test = [traj[8] for traj in traj_central]

# Create orbital plot
p1 = plot(title="Orbital Motion in Central Gravitational Field", 
          xlabel="x (AU)", ylabel="y (AU)", 
          aspect_ratio=:equal, legend=:topright,
          size=(800, 600), dpi=300)

# Orbital trajectories
plot!(p1, x_main, y_main, 
      label="Main Particle", color=:blue, linewidth=2)
plot!(p1, x_test, y_test, 
      label="Test Particle", color=:red, linewidth=2)

# Mark central mass and initial positions
scatter!(p1, [0], [0], label="Central Mass", color=:black, 
            markersize=8, markershape=:star)
scatter!(p1, [q_main[1]], [q_main[2]], label="Initial Positions", 
            color=:green, markersize=6)
scatter!(p1, [q_test[1]], [q_test[2]], label="", color=:green, 
            markersize=6)

# Display and save plot
display(p1)
savefig(p1, "solveT_orbital_motion_test.pdf")
println("  Plot saved as 'solveT_orbital_motion_test.pdf'")