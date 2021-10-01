using Test

using CSV, LinearAlgebra
using ForwardDiff, OrdinaryDiffEq

include("../source/pomin-types.jl")
include("testfunctions.jl")

@testset "All tests:" begin

    @time @testset "idxer tests:" begin 
        include("idxer_test.jl") 
    end

    @time @testset "Integrator tests:" begin 
        include("int_test.jl") 
    end

    @time @testset "Hamiltonian tests:" begin 
        include("hamiltonian_test.jl") 
    end

    @time @testset "Gravitational waveform tests:" begin 
        include("gwsc_test.jl") 
    end

end

nothing