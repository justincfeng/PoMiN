#!/usr/bin/env julia
#-----------------------------------------------------------------------
#
#   Simple PoMiN Test
#   
#   Basic functionality test for the modernized PoMiN codebase
#
#-----------------------------------------------------------------------

# Add the PoMiN source to the path
push!(LOAD_PATH, joinpath(@__DIR__, "..", "src"))

println("🧪 Simple PoMiN Test")
println("===================")

# Test imports
println("\n📦 Testing imports...")
try
    include("../src/pomin.jl")
    println("   ✅ PoMiN module loaded successfully")
catch e
    println("   ❌ Failed to load PoMiN: $e")
    exit(1)
end

# Test type creation
println("\n🔧 Testing type system...")
try
    # Test Parameters creation
    params_rk4 = ParametersRK4((0.0, 10.0), δ=0.1)
    println("   ✅ ParametersRK4 created: $(typeof(params_rk4))")
    
    params_julia = ParametersJulia((0.0, 10.0))
    println("   ✅ ParametersJulia created: $(typeof(params_julia))")
    
    # Test Particles creation
    system = setup_binary_system(1.0, 1.0, 5.0, 0.05)
    println("   ✅ Binary system created: $(typeof(system))")
    println("   • Number of particles: $(length(system.m))")
    println("   • Masses: $(system.m)")
    
catch e
    println("   ❌ Type system test failed: $e")
    exit(1)
end

# Test basic integration
println("\n🚀 Testing integration...")
try
    # Short integration test
    params = ParametersRK4((0.0, 1.0), δ=0.1)
    system = setup_binary_system(1.0, 1.0, 10.0, 0.01)
    
    println("   • Running short RK4 integration...")
    @time solution = solve(system, params)
    
    println("   ✅ Integration completed successfully!")
    println("   • Solution type: $(typeof(solution))")
    println("   • Time steps: $(length(solution.t))")
    println("   • Final time: $(solution.t[end])")
    
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
println("🎉 PoMiN is working correctly!")
