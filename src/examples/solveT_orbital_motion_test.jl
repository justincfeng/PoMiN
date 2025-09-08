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

# Central mass (at origin)
M_central = tpfl(1.0)  # Solar mass units

# Particle masses (equal)
m_particle = tpfl(1e-6)  # Small test masses

# Initial positions (opposite sides of central mass)
r0 = tpfl(5e8)  # Distance from center in AU
q_main = [r0, tpfl(0.0), tpfl(0.0)]      # Main particle at +x
q_test = [-r0, tpfl(0.0), tpfl(0.0)]     # Test particle at -x

# Circular orbital velocities
v_circular = sqrt(M_central / r0)  # Circular velocity
p_main = [tpfl(0.0), m_particle * v_circular, tpfl(0.0)]   # +y momentum
p_test = [tpfl(0.0), -m_particle * v_circular, tpfl(0.0)]  # -y momentum

# Main particle system (interacts gravitationally)
PintMain = pomin.setup_single_particle(m_particle, q_main, p_main, tpfl)

# Test particle system (feels forces but doesn't exert them)
PtestTest = pomin.setup_single_particle(m_particle, q_test, p_test, tpfl)

# Central gravitational potential at origin
Φ_central = ΦCentralMass([tpfl(0.0), tpfl(0.0), tpfl(0.0)], tpfl, M_central)

# One orbital period for comparison
T_orbit = 2*π*sqrt(r0^3/M_central)
tspan = (tpfl(0.0), T_orbit)
tols = tpfl(1e-12)

params = pomin.ParametersJulia(tspan, integrator="Vern9", atol=tols, rtol=tols)

sol_central = pomin.solveT(PintMain, params, PtestTest, Φ_central)

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
scatter!(p1, [0], [0], label="Central Mass", color=:black, markersize=8, markershape=:star)
scatter!(p1, [q_main[1]], [q_main[2]], label="Initial Positions", color=:green, markersize=6)
scatter!(p1, [q_test[1]], [q_test[2]], label="", color=:green, markersize=6)

# Energy conservation check
E_initial = 0.5 * m_particle * v_circular^2 - M_central * m_particle / r0
v_main_final = norm(Z_central[4:6]) / m_particle
r_main_final = norm(q_main_final)
E_final = 0.5 * m_particle * v_main_final^2 - M_central * m_particle / r_main_final
energy_error = abs(E_final - E_initial)/abs(E_initial) * 100

println("  Energy conservation: $(energy_error)%")

# Display and save plot
display(p1)
savefig(p1, "solveT_orbital_motion_test.pdf")
println("  Plot saved as 'solveT_orbital_motion_test.pdf'")