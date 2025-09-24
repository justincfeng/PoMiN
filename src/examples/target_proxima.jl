using LinearAlgebra
using DoubleFloats
tpfl = Double64

# Include relativistic tools
include("../core/physics/Hamiltonians/RelTools.jl")



#-----------------------------------------------------------------------
#
#   FUNCTIONS FOR TARGETING PROBLEM
#
#-----------------------------------------------------------------------

#-----------------------------------------------------------------------
function γV2v(γV)
    #   γ  = sqrt(1+γv^2)
    return γV ./ sqrt(one(typeof(γV[1]))+dot(γV,γV))
end #-------------------------------------------------------------------

#-----------------------------------------------------------------------
function vunit(v) # Returns unit vector in direction of v
    return v/norm(v)
end #-------------------------------------------------------------------

#-----------------------------------------------------------------------
function uvproj(u,v) # Projects vector v onto vector u
    U = vunit(u)
    return v .- dot(U,v) .* U
end #-------------------------------------------------------------------

#-----------------------------------------------------------------------
function basisconstructor(vp,vs) # Constructs basis from two vectors
    epar = vunit(vp)
    e1   = vunit(uvproj(vp,vs))
    e2   = vunit(cross(e1,epar))
    return (epar,e1,e2)
end #-------------------------------------------------------------------

#-----------------------------------------------------------------------
function bconstructor(vecs,b,Θ) # Returns displaced target position
    Xpx0,Xst0,Vpx = vecs

    tpfl = typeof(b)

    ONE,TWO  = one(tpfl),tpfl(2)

    ΔXp0  = Xst0 .- Xpx0

    (epar,e1,e2) = basisconstructor(-ΔXp0,Vpx)

    n = cos(Θ)*e1 + sin(Θ)*e2

    φ = acos(b/norm(ΔXp0))

    return b .* (cos(φ) .* epar + sin(φ) .* n)
end #-------------------------------------------------------------------

#-----------------------------------------------------------------------
function targprox(vecs,vst,bv) # Returns initial conditions for Starchip, Proxima, and q_target_pos_SPJc
    Xpx0,Xst0,Vpx = vecs

    tpfl = typeof(vst)

    ONE,TWO  = one(tpfl),tpfl(2)

    ΔX0  = Xst0 .- (Xpx0 + bv)
    b    = norm(bv)
    rpx  = norm(ΔX0)
    vb   = norm(Vpx)
    
    # From LaTeX: cos(ψ) = -dot(Vpx,ΔX0)/(vb*rpx)
    cosψ = -dot(Vpx,ΔX0) / (vb*rpx)
    sinψ = sqrt(ONE - cosψ^2)
    
    # From LaTeX: sin(ϑ) = (vb/vst) * sin(ψ)
    sinϑ = (vb/vst) * sinψ
    cosϑ = sqrt(ONE - sinϑ^2)
    
    # Construct basis with ΔX0 (not -ΔX0)
    epar = vunit(ΔX0)
    e1   = vunit(uvproj(ΔX0, Vpx))
    
    # Spacecraft velocity: pointing toward target with correct magnitude
    Vst  = vst * (-cosϑ .* epar .+ sinϑ .* e1)
    
    # Velocity difference
    ΔV   = Vst .- Vpx
    
    # Time to closest approach from LaTeX constraint
    tcl  = -dot(ΔX0, ΔV) / dot(ΔV, ΔV)

    return (Xpx0,Xpx0 + bv,Xst0,Vpx,Vst,tcl)
end #-------------------------------------------------------------------

#-----------------------------------------------------------------------
function parfuncs(ics) # Returns functions for Starchip, Proxima, displaced Proxima, and time to closest approach
    Xpx0,Xpxb,Xst0,Vpx,Vst,tcl = ics
    return ( t->(Xst0 .+ (Vst .* t)) , t->(Xpx0 .+ (Vpx .* t)) , 
             t->(Xpxb .+ (Vpx .* t)) , tcl )
end #-------------------------------------------------------------------

#-----------------------------------------------------------------------
#
#   PARAMETER DEFINITIONS
#
#-----------------------------------------------------------------------

#-----------------------------------------------------------------------
#   CONSTANTS
#-----------------------------------------------------------------------
cMKS   = tpfl(299792458.0)
AUMKS  = tpfl(149597870700.0)
lyMKS  = tpfl(9460730472580800.0)
Msol2m = tpfl(1476.67)

c      = one(tpfl)
AU     = AUMKS/Msol2m
ly     = lyMKS/Msol2m 

#-----------------------------------------------------------------------
#   FLAT SPACE TARGETING PROBLEM
#-----------------------------------------------------------------------

# Proxima Centauri
Xpx0 = tpfl.([-9.90338347925777E+12, -7.58520275447461E+12, -2.41494965557852E+13])
γVpx = tpfl.([-3.34779933990201E-5, 7.24198145771899E-5, 6.82553747232694E-5])
Vpx  = γV2v(γVpx)

# Starchip
Xst0 = tpfl.([-2.235888445902370E+07, 8.680664864663420E+07, 2.988257878660500E+07]) # 0.1 AU from Earth, in direction of Proxima's init pos
vst    = tpfl(0.2)*c         # Starchip velocity magnitude

# Target
b0     = tpfl(0.05)*AU               # Target distance
bv     = bconstructor((Xpx0,Xst0,Vpx),b0,tpfl(pi)/tpfl(3)) # Target displacement from Proxima

ics    = targprox( (Xpx0,Xst0,Vpx) , vst , bv )
pfs    = parfuncs(ics)

XstF,XpxF,XbF,tcl = pfs

# pfs[1] is the Starchip position as a function of time
# pfs[2] is the Proxima position as a function of time
# pfs[3] is the displaced target position as a function of time
# pfs[4] is the time to closest approach 

Xpx0,Xpxb,Xst0,Vpx,Vst,tcl = ics

# Store final target position
target_position_final = XbF(tcl)

# Store total travel distance
ΔXst = norm(XstF(tcl) - Xst0)
