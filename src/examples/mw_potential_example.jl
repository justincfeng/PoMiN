#-----------------------------------------------------------------------

using DoubleFloats
using LinearAlgebra
using Printf, Plots
using Statistics

# Include the necessary PoMiN modules
include("../pomin.jl")
include("../core/pomin-types.jl")  # Load type definitions first
include("../core/physics/Hamiltonians/HamTools.jl")  # For Z2q function
include("../core/physics/external_potentials/external.jl")

tpfl = Double64

# Import the gradient operator from pomin module
∂ = pomin.∂

#-----------------------------------------------------------------------

# Solar offset values (in solar mass units):
origin_x = tpfl(-1.708859462494220e17)  # Solar x-offset
origin_z = tpfl(4.346342845091530e14)   # Solar z-offset
origin_y = tpfl(0.0)                    # Solar y-offset

xo = [origin_x, origin_y, origin_z]

# Velocity conversion factor
km_s = tpfl(1000.0 / 299792458.0)

# Initial velocity of the Sun
vpec    = km_s .* tpfl.([11.1, 12.24, 7.25])  # Solar peculiar motion
v_LSR   = tpfl.([0.0, 220.0 * km_s, 0.0])     # LSR circular motion
v_total = v_LSR + vpec                         # Total velocity

γ       = one(tpfl)/sqrt(one(tpfl) - dot(v_total, v_total))

u_total = γ .* v_total

# Milky Way potential 
Φ = ΦMilkyWay(tpfl, xo)

T_orbit = tpfl(2*π) * norm(xo) / norm(v_total)

# Solar mass only
M_sun = one(tpfl)   # Solar mass in solar mass units
m     = [M_sun]     # Just the Sun

# Construct the potential energy function
V     = UConstructor(Φ, m, zeros(tpfl,3))

# Construct initial data
q0    = [zeros(tpfl,3)]
p0    = [m[1] .* u_total] 

# Construct the particle system (just the sun)
system = pomin.Particles(m, q0, p0)

tspan = (0.0, T_orbit * 2.0)  # Simulate for 2 complete orbits
params = pomin.ParametersJulia(tspan)

sol = pomin.solve(system, params, V)

#-----------------------------------------------------------------------

# Extract trajectory functions
qs = t->sol(t)[1:3]
vs = t->sol(t)[4:6] ./ m[1]

# Unit conversions
kpc = 3.0857e19 / 1476.67  # 1 kpc in solar mass units

# Extract data
times = sol.t
positions = [qs(t) for t in times]
x_traj = [pos[1] for pos in positions]
y_traj = [pos[2] for pos in positions]

# Energy conservation
total_energy = [0.5 * dot(vs(t), vs(t)) + Φ(qs(t)) for t in times]
energy_variation = (total_energy .- total_energy[1]) ./ 
                    abs(total_energy[1])

println("Energy conservation: max variation = 
        $(round(maximum(abs.(energy_variation)) * 100, digits=6))%")

println("Average speed of Sun = 
        $(round(mean([norm(vs(t)) for t in times]) * 100, digits=6))%")

# Plots
p1 = plot(x_traj ./ kpc, y_traj ./ kpc, 
          label="Solar Orbit", linewidth=2, color=:orange,
          title="Solar Orbit", xlabel="x (kpc)", ylabel="y (kpc)",
          aspect_ratio=:equal)
scatter!(p1, [x_traj[1] / kpc], [y_traj[1] / kpc], 
         label="Start", markersize=4, color=:green)

p2 = plot(times ./ T_orbit, energy_variation .* 100,
          label="Energy Variation", linewidth=2, color=:green,
          title="Energy Conservation", xlabel="Periods", 
          ylabel="ΔE/E₀ (%)")
hline!(p2, [0.0], label="Conserved", linestyle=:dash, color=:red)

combined_plot = plot(p1, p2, layout=(1,2), size=(800,400))
display(combined_plot)
savefig(combined_plot, "mw_sun_orbit.pdf")

#-----------------------------------------------------------------------
