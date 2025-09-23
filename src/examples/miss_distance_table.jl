using LinearAlgebra
using Printf

# Include the target position calculation
include("target_proxima.jl")

# Include the main PoMiN module
include("../pomin.jl")
using .pomin

println("="^80)
println("MISS DISTANCE ANALYSIS TABLE")
println("="^80)

# Get the target position from the flat space script
XtarF = target_position_final

# Initialize storage for results
bodies = ["Sun", "Alpha Centauri", "Jupiter", "Earth", "Proxima", "Moon", "Mars"]
newtonian_misses = Float64[]
relativistic_misses = Float64[]

# Common parameters (from the original script)
tpfl = Double64
tcl = tpfl(1.35858929869399063971729181215165883e+14)  # Time to closest approach

# Spacecraft initial conditions
qchip = tpfl.([-2.235888445902370E+07, 8.680664864663420E+07, 2.988257878660500E+07])
pchip_rel = tpfl.([1.0e-10 * 0.2 * 1.0, 1.0e-10 * 0.2 * 1.0, 1.0e-10 * (-0.2) * 1.0])
pchip_newt = pchip_rel

# Create test particles
Pchip_rel0 = pomin.Particles([1.0e-10], [qchip], [pchip_rel])
Pchip_newt0 = pomin.Particles([1.0e-10], [qchip], [pchip_newt])

# Integration parameters
params = (tspan=(tpfl(0.0), tcl), reltol=1e-15, abstol=1e-17)

println("Calculating miss distances for each celestial body...")
println()

# Function to calculate miss distance for a given particle system
function calculate_miss_distance(particle_system, body_name)
    try
        # Newtonian case
        solN = pomin.solve(particle_system, params; testparticles=Pchip_newt0, Newtonian=true)
        zendN = solN(tcl)
        qmissN = zendN[7:9] - XtarF
        dmissN = norm(qmissN)
        
        # Relativistic case  
        sol = pomin.solve(particle_system, params; testparticles=Pchip_rel0)
        zend = sol(tcl)
        qmiss = zend[7:9] - XtarF
        dmiss = norm(qmiss)
        
        println("$body_name:")
        println("  Newtonian miss distance:    $(dmissN)")
        println("  Relativistic miss distance: $(dmiss)")
        println("  Ratio (N/R):               $(dmissN/dmiss)")
        println()
        
        return dmissN, dmiss
    catch e
        println("Error calculating $body_name: $e")
        return NaN, NaN
    end
end

# Note: You would need to define the particle systems (Psol, Palpha, etc.) 
# from the original script to make this work. For now, this shows the structure.

println("Table Structure:")
println(@sprintf("%-15s | %-18s | %-19s | %-10s", "Body", "Newtonian Miss", "Relativistic Miss", "Ratio"))
println("-"^70)

# This would be populated with actual calculations:
# for (i, body) in enumerate(bodies)
#     dmissN, dmiss = calculate_miss_distance(particle_systems[i], body)
#     println(@sprintf("%-15s | %-18.6e | %-19.6e | %-10.3f", body, dmissN, dmiss, dmissN/dmiss))
# end

println("Note: To complete this table, you need to:")
println("1. Define all particle systems (Psol, Palpha, Pjup, PEarth, PProx, PMoon, PMars)")
println("2. Run the calculate_miss_distance function for each system")
println("3. The results will be automatically formatted in the table above")
