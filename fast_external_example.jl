using DoubleFloats
using LinearAlgebra
using Printf

# Include the necessary PoMiN modules
include("src/pomin.jl")
include("src/core/pomin-types.jl")
include("src/core/physics/Hamiltonians/HamTools.jl")
include("src/core/physics/external_potentials/external.jl")

println("🚀 Fast External Potential Example")
println("=" ^ 50)

# Use circular orbit with reduced velocity for relativistic case
r0 = 10.0
M = 1.0
velocity_factor = 0.2
ve = velocity_factor * sqrt(M/r0)
Te = 2*π*sqrt(r0^3/M)

V = ΦCentralMass(zeros(3), Double64, 1.0)
m = [1e-10, 1e-10]

q0 = [[r0, 0.0, 0.0], [-r0, 0.0, 0.0]]
p0 = [[0.0, m[1]*ve, 0.0], [0.0, -m[2]*ve, 0.0]]

system = pomin.Particles(m, q0, p0)

# Use faster integration parameters
tspan = (0.0, Te/4)  # Shorter time span for testing
params = pomin.ParametersJulia(tspan, atol=1e-6, rtol=1e-6)  # Relaxed tolerances

println("Starting integration...")
@time sol = pomin.solve(system, params, V)

# Quick statistics without plotting
n_steps = length(sol.t)
x1_traj = [sol.u[i][1] for i in 1:n_steps]
y1_traj = [sol.u[i][2] for i in 1:n_steps]
x2_traj = [sol.u[i][4] for i in 1:n_steps]
y2_traj = [sol.u[i][5] for i in 1:n_steps]

separation = [sqrt((x1_traj[i] - x2_traj[i])^2 + (y1_traj[i] - y2_traj[i])^2) for i in 1:n_steps]

println("\n=== Quick Results ===")
println("Integration time: $(sol.t[end]) time units")
println("Number of time steps: $n_steps")
println("Initial separation: $(separation[1])")
println("Final separation: $(separation[end])")
println("Max separation: $(maximum(separation))")
println("Min separation: $(minimum(separation))")

println("\n✅ Fast example completed!")
