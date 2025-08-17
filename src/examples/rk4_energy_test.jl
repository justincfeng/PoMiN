#!/usr/bin/env julia

"""
RK4 Energy Conservation Diagnostic
Test if the RK4 integrator is conserving energy properly.
"""

using Printf
using DoubleFloats
include("../pomin.jl")

function test_energy_conservation()
    println("=== RK4 Energy Conservation Test ===")
    
    # Simple test case
    m1, m2 = 1.0, 1.0
    p = 10.0
    b = 10.0
    dx = 100.0
    
    # Set up system
    system = pomin.setup_scattering(m1, m2, p, b, dx, tpfl=Double64)
    
    # Short integration time for diagnostic
    t_flight = Double64(10.0)
    δ_initial = Double64(0.01)
    courant = Double64(0.001)
    
    # Solve with many saved points to check energy
    sol = pomin.solve(system, pomin.ParametersRK4((Double64(0.0), t_flight), 
                     δ=δ_initial, courant=courant, Nrec=100))
    
    println("Number of time points saved: $(length(sol.t))")
    
    # Check energy conservation
    initial_energy = pomin.HamPM.H(sol.z[1], [m1, m2])
    final_energy = pomin.HamPM.H(sol.z[end], [m1, m2])
    
    println("Initial energy: $(initial_energy)")
    println("Final energy:   $(final_energy)")
    println("Energy change:  $(final_energy - initial_energy)")
    println("Relative change: $((final_energy - initial_energy) / initial_energy * 100)%")
    
    # Check energy at several points
    println("\nEnergy evolution:")
    for i in 1:min(10, length(sol.t))
        E = pomin.HamPM.H(sol.z[i], [m1, m2])
        ΔE = (E - initial_energy) / initial_energy * 100
        println("  t=$(round(sol.t[i], digits=3)), E=$(round(E, digits=6)), ΔE=$(round(ΔE, digits=8))%")
    end
    
    # Check if particles are moving in expected direction
    println("\nParticle motion check:")
    q1_initial = [sol.z[1][1], sol.z[1][2], sol.z[1][3]]
    q1_final = [sol.z[end][1], sol.z[end][2], sol.z[end][3]]
    p1_initial = [sol.z[1][7], sol.z[1][8], sol.z[1][9]]
    p1_final = [sol.z[end][7], sol.z[end][8], sol.z[end][9]]
    
    println("  Particle 1 initial position: $(q1_initial)")
    println("  Particle 1 final position:   $(q1_final)")
    println("  Particle 1 initial momentum: $(p1_initial)")
    println("  Particle 1 final momentum:   $(p1_final)")
    
    # Check if particle moved in +x direction (as expected)
    Δx = q1_final[1] - q1_initial[1]
    println("  Particle 1 x-displacement: $(Δx)")
    if Δx > 0
        println("  ✅ Particle 1 moved in +x direction (expected)")
    else
        println("  ❌ Particle 1 moved in -x direction (WRONG!)")
    end
    
    return sol
end

# Run the test
test_energy_conservation()
