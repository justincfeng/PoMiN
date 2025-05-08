using ForwardDiff
using LinearAlgebra
using DoubleFloats
using Test

include("pomin-types.jl")
include("idxer.jl")
include("HamPM.jl")

# Test setup for massless particles
function test_massless_hamiltonian(n::Int=10)
    # Set up system parameters
    d = 3  # spatial dimensions
    tol = 1e-12  # tolerance for comparisons
    
    println("Running $n test cases...\n")
    
    # Helper function for H3F0
    function ybaf0(mb::Real, qb::RealVec, qa::RealVec, pb::RealVec, pa::RealVec)
        Θba = Θabf(qb,qa,pb)
        return sqrt(mb^2+Θba^2)/Enf(mb,psf(pb))
    end
    
    # Helper function for original H3F
    function H3F(d::Int, m::RealVec, Z::RealVec)
        tpfl = typeof(Z[1])
        n = length(m)
        qa, qb, pa, pb = [zeros(tpfl,d) for _ = 1:4]
        psa, psb, Ena, Enb = [zero(tpfl) for _ = 1:4]
        rab, yba, Θab, Θba, Ξab = [zero(tpfl) for _ = 1:5]
        H3 = zero(tpfl)
        
        for a=1:n
            qa = Z2q(n,d,a,Z)
            pa = Z2p(n,d,a,Z)
            psa = psf(pa)
            Ena = Enf(m[a],psa)
            
            for b=1:n
                if b!=a
                    qb = Z2q(n,d,b,Z)
                    pb = Z2p(n,d,b,Z)
                    psb = psf(pb)
                    Enb = Enf(m[b],psb)
                    rab = rf(qa,qb)
                    yba = ybaf(m[b],qb,qa,pb,pa)
                    Θab = Θabf(qa,qb,pa)
                    Θba = Θabf(qb,qa,pb)
                    Ξab = Ξabf(pa,pb)
                    
                    H3 += ( 2*(-2*(Ξab*Θba)^2 - 2*Θab*Θba*Ξab*psb - (Θab*psb)^2 + psb*Ξab^2)/Enb^2 + 
                            2*(psa*Θba^2 - (Θab*Θba)^2 + 2*Θab*Θba*Ξab - Ξab^2 + psb*Θab^2) +
                            yba*(3*psa*Θba^2 - (Θab*Θba)^2 + 8*Θab*Θba*Ξab - psa*psb + 3*psb*Θab^2)
                            ) / (4*Ena*Enb*rab*yba*(yba+1)^2)
                end
            end
        end
        return H3
    end
    
    # Helper function for new H3F0
    function H3F0(d::Int, m::RealVec, Z::RealVec)
        tpfl = typeof(Z[1])
        n = length(m)
        qa, qb, pa, pb = [zeros(tpfl,d) for _ = 1:4]
        psa, psb, Ena, Enb = [zero(tpfl) for _ = 1:4]
        rab, yba, Θab, Θba, Ξab = [zero(tpfl) for _ = 1:5]
        H3 = zero(tpfl)
        
        for a=1:n
            qa = Z2q(n,d,a,Z)
            pa = Z2p(n,d,a,Z)
            psa = psf(pa)
            Ena = Enf(m[a],psa)
            
            for b=1:n
                if b!=a
                    qb = Z2q(n,d,b,Z)
                    pb = Z2p(n,d,b,Z)
                    psb = psf(pb)
                    Enb = Enf(m[b],psb)
                    rab = rf(qa,qb)
                    yba = ybaf0(m[b],qb,qa,pb,pa)
                    Θab = Θabf(qa,qb,pa)
                    Θba = Θabf(qb,qa,pb)
                    Ξab = Ξabf(pa,pb)
                    
                    H3 += (3*Enb*psb*Θab^2 - Enb*psa*(psb - 3*Θba^2) + Enb*Θab*Θba*(-(Θab*Θba) + 8*Ξab) + 
                           2*psa*psb*abs(Θba) - 2*psb*Θab^2*abs(Θba) - 4*Ξab^2*abs(Θba))/(4*Ena*rab*(Enb + abs(Θba))^2)
                end
            end
        end
        return H3
    end
    
    @testset "Massless Particle Tests" begin
        for i in 1:n
            # Both particles are massless
            m = RealVec{Double64}([0.0, 0.0])  # Two massless particles
            
            # Initial positions with finite separation
            q1 = RealVec{Double64}(zeros(d))  # First particle at origin
            q2 = RealVec{Double64}(rand(d))   # Random position for second particle
            
            # Initial momenta (non-orthogonal)
            p1 = RealVec{Double64}(rand(d))  # Random momentum for first particle
            p2 = RealVec{Double64}(rand(d))  # Random momentum for second particle
            
            # Combine into phase space vector Z
            Z = RealVec{Double64}(vcat(vcat(q1, q2), vcat(p1, p2)))
            
            # Compute Hamiltonians
            h_old = H3F(d, m, Z)
            h_new = H3F0(d, m, Z)
            
            # Compute gradients
            grad_old = Float64.(ForwardDiff.gradient(z -> H3F(d, m, z), Z))
            grad_new = Float64.(ForwardDiff.gradient(z -> H3F0(d, m, z), Z))
            
            # Debug info
            println("Case $i:")
            println("  H3F  = $(Float64(h_old))")
            println("  H3F0 = $(Float64(h_new))")
            
            # Test Hamiltonian values match
            @test isapprox(Float64(h_old), Float64(h_new), rtol=tol)
            
            # Test gradients match
            @test all(isapprox.(grad_old, grad_new, rtol=tol))
        end
    end
end

println("Testing Hamiltonian for massless particles...")
test_massless_hamiltonian(100)