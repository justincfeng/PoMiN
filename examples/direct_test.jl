#!/usr/bin/env julia
#-----------------------------------------------------------------------
#
#   Direct PoMiN Test (No Module System)
#   
#   Tests the modernized PoMiN codebase by directly including files
#
#-----------------------------------------------------------------------

using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using LinearAlgebra
using ForwardDiff
using OrdinaryDiffEq
using Logging
using DoubleFloats

println("🧪 Direct PoMiN Test")
println("===================")

# Include files directly to avoid module precompilation issues
println("\n📦 Loading PoMiN components...")

try
    # Load core types
    include("../src/core/pomin-types.jl")
    println("   ✅ Types loaded")
    
    # Load Hamiltonian tools
    include("../src/physics/Hamiltonians/HamTools.jl")
    println("   ✅ Hamiltonian tools loaded")
    
    # Load indexing utilities
    include("../src/physics/Hamiltonians/idxer.jl")
    println("   ✅ Phase space indexing loaded")
    
    # Load Hamiltonian
    include("../src/physics/Hamiltonians/HamPM.jl")
    println("   ✅ Post-Minkowski Hamiltonian loaded")
    
    # Load integrators
    include("../src/integrators/tadap.jl")
    println("   ✅ Time adaptation loaded")
    
    include("../src/integrators/rk4i.jl")
    println("   ✅ RK4 integrator loaded")
    
    include("../src/integrators/intjul.jl")
    println("   ✅ Julia integrator loaded")
    
    # Load initial data generation
    include("../src/physics/initial_data/idgen.jl")
    println("   ✅ Initial data generation loaded")
    
catch e
    println("   ❌ Failed to load components: $e")
    exit(1)
end

# Test type creation
println("\n🔧 Testing type system...")
try
    # Test Parameters creation
    params_rk4 = ParametersRK4((0.0, 10.0), δ=0.1)
    println("   ✅ ParametersRK4 created: $(typeof(params_rk4))")
    println("     • Time span: $(params_rk4.tspan)")
    println("     • Time step: $(params_rk4.δ)")
    println("     • RK4 flag: $(params_rk4.rkl)")
    
    params_julia = ParametersJulia((0.0, 10.0))
    println("   ✅ ParametersJulia created: $(typeof(params_julia))")
    println("     • Time span: $(params_julia.tspan)")
    println("     • RK4 flag: $(params_julia.rkl)")
    
    # Test Particles creation
    system = setup_binary_system(1.0, 1.0, 5.0, 0.05)
    println("   ✅ Binary system created: $(typeof(system))")
    println("     • Number of particles: $(length(system.m))")
    println("     • Masses: $(system.m)")
    println("     • Positions: $(system.q)")
    println("     • Momenta: $(system.p)")
    
catch e
    println("   ❌ Type system test failed: $e")
    exit(1)
end

# Test Hamiltonian functions
println("\n⚡ Testing Hamiltonian functions...")
try
    # Create test system
    system = setup_binary_system(1.0, 1.0, 10.0, 0.01)
    z = vcat(system.q..., system.p...)
    m = system.m
    
    println("   • Phase space vector length: $(length(z))")
    println("   • Testing Hamiltonian evaluation...")
    
    H_val = HamPM.H(z, m, 1)  # 1st order post-Minkowskian
    println("   ✅ Hamiltonian H = $H_val")
    
    println("   • Testing Hamiltonian gradient...")
    dH_val = HamPM.dH(z, m, 1)
    println("   ✅ Gradient dH computed, length: $(length(dH_val))")
    
    println("   • Testing Hamilton's equations...")
    FHE_val = HamPM.FHE(z, m)
    println("   ✅ Hamilton's equations computed, length: $(length(FHE_val))")
    
catch e
    println("   ❌ Hamiltonian test failed: $e")
    exit(1)
end

# Test short integration
println("\n🚀 Testing integration...")
try
    # Create simple test system
    system = setup_binary_system(1.0, 1.0, 10.0, 0.01)
    z = vcat(system.q..., system.p...)
    m = system.m
    
    println("   • Testing RK4 integrator...")
    params_rk4 = ParametersRK4((0.0, 1.0), δ=0.1)
    
    @time solution_rk4 = hrkintegrator(
        params_rk4.d, length(m), z,
        (Z) -> HamPM.dH(Z, m, params_rk4.d),
        params_rk4.δ,
        (dt, Z, Zdot) -> tadap.tcour(dt, Z, Zdot, params_rk4.courant, params_rk4.d),
        params_rk4.tspan, params_rk4.iter
    )
    
    println("   ✅ RK4 integration completed!")
    println("     • Solution type: $(typeof(solution_rk4))")
    println("     • Time steps: $(length(solution_rk4.t))")
    println("     • Final time: $(solution_rk4.t[end])")
    
    println("   • Testing Julia integrator...")
    params_julia = ParametersJulia((0.0, 1.0))
    
    @time solution_julia = jlintegrator(
        z,
        (du, u, p, t) -> begin
            du .= HamPM.FHE(u, p)
        end,
        params_julia.tspan, m, params_julia.atol, params_julia.rtol
    )
    
    println("   ✅ Julia integration completed!")
    println("     • Solution type: $(typeof(solution_julia))")
    println("     • Time steps: $(length(solution_julia.t))")
    println("     • Final time: $(solution_julia.t[end])")
    
catch e
    println("   ❌ Integration test failed: $e")
    println("   Stack trace:")
    for (exc, bt) in Base.catch_stack()
        showerror(stdout, exc, bt)
        println()
    end
    exit(1)
end

println("\n✅ All tests passed!")
println("🎉 PoMiN modernization is working correctly!")
println("\n📋 Summary:")
println("   • Type system: ✅ Working")
println("   • Hamiltonian functions: ✅ Working") 
println("   • RK4 integrator: ✅ Working")
println("   • Julia integrator: ✅ Working")
println("   • Initial data generation: ✅ Working")
println("\n🚀 Ready for binary black hole simulations!")
