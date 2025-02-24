using LinearAlgebra
using Test

# Include type definitions
include("../src/pomin-types.jl")

# Include the source file
include("../src/HamPM.jl")

# Test massless particle scattering
function test_massless_scattering()
    # Set up initial conditions for two massless particles
    n = 2  # number of particles
    d = 3  # dimensions
    m = zeros(Float64, n)  # massless particles
    
    # Initial positions: particles start separated in x direction
    q1 = [1.0, 0.0, 0.0]
    q2 = [-1.0, 0.0, 0.0]
    
    # Initial momenta: particles moving towards each other
    p1 = [-1.0, 0.0, 0.0]
    p2 = [1.0, 0.0, 0.0]
    
    # Combine into state vector
    Z = vcat(q1, q2, p1, p2)
    
    # Calculate initial Hamiltonian
    H_init = H3m0(d, m, Z)
    
    # Calculate derivatives
    dH = dH3m0(d, m, Z)
    
    # Test energy conservation
    @test abs(sum(dH .* Z)) < 1e-10
    
    # Test momentum conservation
    for i in 1:d
        total_momentum = sum(Z[2*n*d-d+i:d:2*n*d])
        @test abs(total_momentum) < 1e-10
    end
    
    # Test that derivatives are non-zero (interaction is happening)
    @test norm(dH) > 0
    
    println("Initial Hamiltonian: ", H_init)
    println("Derivative norm: ", norm(dH))
end

# Run tests
test_massless_scattering()
