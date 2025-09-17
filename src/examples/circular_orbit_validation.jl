using LinearAlgebra
using Plots
using Statistics
using DoubleFloats

include("../pomin.jl")

#-----------------------------------------------------------------------
#   KEPLERIAN ORBIT FUNCTION
#-----------------------------------------------------------------------
"""
    keplerian_orbit(a,M,e,periapsis=false)

Returns the initial conditions for a Keplerian orbit at apoapsis 
(maximum separation).
"""
function keplerian_orbit(a,M,e,periapsis=false)
    if periapsis
        r0 = a*(1-e)  # Periapsis
    else
        r0 = a*(1+e)  # Apoapsis
    end
    ve = sqrt(abs( M*(2.0/r0-1.0/a) ))
    l  = ve*r0            # Angular momentum
    ϵ = ve^2/2 - M/r0     # Total energy
    Te = 2*π*sqrt(a^3/M)
    return (r0,ve,l,ϵ,Te)
end

#-----------------------------------------------------------------------
#   PARAMETERS
#-----------------------------------------------------------------------

tpfl  =  Double64

ν0 = zero(tpfl)
ν1 = one(tpfl)
ν2 = tpfl(2)

# Circular orbit setup
a0 = tpfl(1e6)    # Orbital radius
e0 = ν0           # Eccentricity
M  = ν1           # Central mass

(r0,ve,l,ϵ,Te) = keplerian_orbit(a0,M,e0,false)

Te = ν2*π*sqrt(a0^3/M)     # Orbital period

#-----------------------------------------------------------------------
#   MASSIVE CENTRAL PARTICLE
#-----------------------------------------------------------------------

qC   = zeros(tpfl,3)
pC   = zeros(tpfl,3)

# Particle object for central mass
PC   = pomin.setup_single_particle(M, qC, pC, tpfl)

#-----------------------------------------------------------------------
#   LIGHT PARTICLE
#-----------------------------------------------------------------------

m    = tpfl(1e-6)    # Light particle mass

γL = ν1/sqrt(ν1-ve^2)

qL   = [a0 , ν0 , ν0 ]
pL   = [ν0 , m*γL*ve , ν0 ]

# Particle object for light particle
PL   = pomin.setup_single_particle(m, qL, pL, tpfl)

PCombined = pomin.merge_particle_systems(PC, PL)

#-----------------------------------------------------------------------
#   PARAMETERS
#-----------------------------------------------------------------------

tspan = (zero(tpfl), Te)

tols  = tpfl(1e-16)

params       = pomin.ParametersJulia( tspan, integrator="Vern9", 
                                      atol=tols, rtol=tols)      

#-----------------------------------------------------------------------
#   SOLVE
#-----------------------------------------------------------------------

solN         = pomin.solve(PCombined, params; Newtonian=true)
sol          = pomin.solve(PCombined, params)
solT         = pomin.solve(PC, params; testparticles=PL)

#-----------------------------------------------------------------------
#   PLOTTING
#-----------------------------------------------------------------------

println("Creating comparative orbital trajectory plots...")

# Function to extract trajectories from different solution types
function extract_trajectories(solution, solver_name)
    n_pts = length(solution.t)
    
    # Handle different solution formats
    if hasfield(typeof(solution), :z)
        # RK4 solutions use .z
        if solver_name == "Test Particle"
            # Test particle: elements 7-9 are spatial vectors for light particle
            xC = [solution.z[i][1] for i in 1:n_pts]
            yC = [solution.z[i][2] for i in 1:n_pts]
            zC = [solution.z[i][3] for i in 1:n_pts]
            xL = [solution.z[i][7] for i in 1:n_pts]
            yL = [solution.z[i][8] for i in 1:n_pts]
            zL = [solution.z[i][9] for i in 1:n_pts]
        else
            # Standard format for other solvers
            xC = [solution.z[i][1] for i in 1:n_pts]
            yC = [solution.z[i][2] for i in 1:n_pts]
            zC = [solution.z[i][3] for i in 1:n_pts]
            xL = [solution.z[i][4] for i in 1:n_pts]
            yL = [solution.z[i][5] for i in 1:n_pts]
            zL = [solution.z[i][6] for i in 1:n_pts]
        end
    else
        # Julia ODE solutions use .u
        if solver_name == "Test Particle"
            # Test particle: elements 7-9 are spatial vectors for light particle
            xC = [solution.u[i][1] for i in 1:n_pts]
            yC = [solution.u[i][2] for i in 1:n_pts]
            zC = [solution.u[i][3] for i in 1:n_pts]
            xL = [solution.u[i][7] for i in 1:n_pts]
            yL = [solution.u[i][8] for i in 1:n_pts]
            zL = [solution.u[i][9] for i in 1:n_pts]
        else
            # Standard format for other solvers
            xC = [solution.u[i][1] for i in 1:n_pts]
            yC = [solution.u[i][2] for i in 1:n_pts]
            zC = [solution.u[i][3] for i in 1:n_pts]
            xL = [solution.u[i][4] for i in 1:n_pts]
            yL = [solution.u[i][5] for i in 1:n_pts]
            zL = [solution.u[i][6] for i in 1:n_pts]
        end
    end
    
    # Calculate separations
    separations = [sqrt((xL[i] - xC[i])^2 + (yL[i] - yC[i])^2 + (zL[i] - zC[i])^2) for i in 1:n_pts]
    
    return (xC, yC, xL, yL, separations, solution.t, n_pts)
end

# Extract trajectories from all three solutions
(xC_N, yC_N, xL_N, yL_N, sep_N, t_N, n_N) = extract_trajectories(solN, "Newtonian")
(xC_PM, yC_PM, xL_PM, yL_PM, sep_PM, t_PM, n_PM) = extract_trajectories(sol, "Post-Minkowskian")
(xC_T, yC_T, xL_T, yL_T, sep_T, t_T, n_T) = extract_trajectories(solT, "Test Particle")

# Create separate orbital trajectory plots for each solver
p1 = plot(title="Newtonian Orbital Motion", xlabel="x", ylabel="y",
          aspect_ratio=:equal, grid=true, legend=:topright)
plot!(p1, xC_N, yC_N, label="Central Mass", linewidth=3, color=:red, marker=:circle, markersize=2)
plot!(p1, xL_N, yL_N, label="Light Particle", linewidth=2, color=:blue)
scatter!(p1, [xL_N[1]], [yL_N[1]], label="Start", markersize=8, color=:green)
scatter!(p1, [xL_N[end]], [yL_N[end]], label="End", markersize=8, color=:orange)

p2 = plot(title="Post-Minkowskian Orbital Motion", xlabel="x", ylabel="y",
          aspect_ratio=:equal, grid=true, legend=:topright)
plot!(p2, xC_PM, yC_PM, label="Central Mass", linewidth=3, color=:red, marker=:circle, markersize=2)
plot!(p2, xL_PM, yL_PM, label="Light Particle", linewidth=2, color=:darkblue)
scatter!(p2, [xL_PM[1]], [yL_PM[1]], label="Start", markersize=8, color=:green)
scatter!(p2, [xL_PM[end]], [yL_PM[end]], label="End", markersize=8, color=:orange)

p3 = plot(title="Test Particle Orbital Motion", xlabel="x", ylabel="y",
          aspect_ratio=:equal, grid=true, legend=:topright)
plot!(p3, xC_T, yC_T, label="Central Mass", linewidth=3, color=:red, marker=:circle, markersize=2)
plot!(p3, xL_T, yL_T, label="Light Particle", linewidth=2, color=:green)
scatter!(p3, [xL_T[1]], [yL_T[1]], label="Start", markersize=8, color=:green)
scatter!(p3, [xL_T[end]], [yL_T[end]], label="End", markersize=8, color=:orange)

# Calculate relative differences from theoretical radius a0
rel_diff_N = [(r - a0)/a0 * 100 for r in sep_N]
rel_diff_PM = [(r - a0)/a0 * 100 for r in sep_PM]
rel_diff_T = [(r - a0)/a0 * 100 for r in sep_T]

p4 = plot(title="Relative Difference from Theoretical Radius a₀", xlabel="Time", ylabel="(r - a₀)/a₀ (%)",
          grid=true, legend=:topright)
plot!(p4, t_N, rel_diff_N, label="Newtonian", linewidth=2, color=:blue)
plot!(p4, t_PM, rel_diff_PM, label="Post-Minkowskian", linewidth=2, color=:red)
plot!(p4, t_T, rel_diff_T, label="Test Particle", linewidth=2, color=:green)

# Combine all plots
final_plot = plot(p1, p2, p3, p4, layout=(2,2), size=(1200, 1000))

# Save plot
savefig(final_plot, "circular_orbit_comparison.pdf")
println("Plot saved as circular_orbit_comparison.pdf")

# Calculate statistics for printing
r_mean_N = mean(sep_N)
r_mean_PM = mean(sep_PM)
r_mean_T = mean(sep_T)


# Print comparative statistics
println("\nComparative Orbital Statistics:")
println("="^50)
println("Newtonian Solver:")
println("  Integration points: $n_N")
println("  Final time: $(t_N[end])")
println("  Mean orbital radius: $r_mean_N")
println("  Max relative difference from a₀: $(maximum(abs.(rel_diff_N)))%")

println("\nPost-Minkowskian Solver:")
println("  Integration points: $n_PM")
println("  Final time: $(t_PM[end])")
println("  Mean orbital radius: $r_mean_PM")
println("  Max relative difference from a₀: $(maximum(abs.(rel_diff_PM)))%")

println("\nTest Particle Solver:")
println("  Integration points: $n_T")
println("  Final time: $(t_T[end])")
println("  Mean orbital radius: $r_mean_T")
println("  Max relative difference from a₀: $(maximum(abs.(rel_diff_T)))%")

println("\nTheoretical Values:")
println("  Orbital radius: $a0")
println("  Orbital period: $Te")
