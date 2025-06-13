#-----------------------------------------------------------------------
#
#   PoMiN: Main solver
#
#-----------------------------------------------------------------------

function pomin_solv(system::Particles, params::Parameters)
    m = system.m
    z = vcat(system.q...,system.p...)
    
    if params.jli[1]
        include("integrators/intjul.jl")
        return jlintegratorfull(z,HamPM.FHE,m,params.tspan)    
    elseif params.rkl[1]
        include("integrators/rk4i.jl")
        include("integrators/tadap.jl")
        return rkl_solv(z,HamPM.FHE,m,params.tspan)
    end
end #-------------------------------------------------------------------
