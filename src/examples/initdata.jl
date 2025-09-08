include("../pomin.jl")
using .pomin
using LinearAlgebra, Printf, Plots
using DoubleFloats
tpfl  = Double64

include("../core/initial_data/idgen.jl")

testInitVectorsForUniformity(0.001,7)