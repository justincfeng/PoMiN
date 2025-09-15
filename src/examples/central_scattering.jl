#-----------------------------------------------------------------------
#
#   3-BODY CLOSEST APPROACH CALCULATION
#   Modernized version compatible with current PoMiN framework
#
#-----------------------------------------------------------------------

#!/usr/bin/env julia

include("../pomin.jl")
using .pomin
using LinearAlgebra, DoubleFloats, Printf
using Plots

include("target_proxima.jl")

# Import Z2q function from HamTools
using .pomin: Z2q
tpfl  =  Double64  # Use Double64 to match target_proxima.jl

#-----------------------------------------------------------------------
#   INITIAL DATA SETUP FOR SPACECRAFT
#-----------------------------------------------------------------------

Dist = norm(Xpx0) ./ tpfl(2)

b    = AU

# Mass and position
mchip = tpfl(1.0E-3)
qchip = [-Dist, b, zero(tpfl) ]

γst = one(tpfl)/sqrt(one(tpfl)-(vst/c)^2)

Vst = [ vst , zero(tpfl), zero(tpfl) ]

# Momentum scaled by same factor as mass to preserve velocity
pchip = (mchip*γst) .* (Vst)

# Particle object for spacecraft
Pchip = pomin.setup_single_particle(mchip, qchip, pchip, tpfl)

#-----------------------------------------------------------------------
#   INITIAL DATA SETUP FOR PROXIMA CENTAURI
#-----------------------------------------------------------------------

# Mass, position, momentum
mProx = tpfl(0.1221)
qProx = zeros(tpfl,3)
pProx = zeros(tpfl,3)

# Particle object for Proxima Centauri
PProx = pomin.setup_single_particle(mProx, qProx, pProx, tpfl)

m2    = tpfl(1.0E-9)

aux_pos = [Dist, 0.0, 0.0]
aux_mom = [0.0, -m2 * tpfl(1.0E-9) , 0.0]

Paux = pomin.add_particle(PProx, m2, aux_pos, aux_mom)

#-----------------------------------------------------------------------
#   
#-----------------------------------------------------------------------

tspan = (zero(tpfl),2*Dist/tpfl(0.2))

δ     = tpfl(7.31E+8) # about one hour
tols  = tpfl(1e-16)

Pint = Paux
Ptest = Pchip

params       = pomin.ParametersJulia( tspan, integrator="Vern9", 
                                     atol=tols, rtol=tols)

sol          = pomin.solveT(Pint, params, Ptest)

sol.u[end]-sol.u[1]

#-----------------------------------------------------------------------
#   PLOTTING ROUTINE
#-----------------------------------------------------------------------

# Debug: Check solution structure
println("Solution type: ", typeof(sol))
println("Solution fields: ", fieldnames(typeof(sol)))
if hasfield(typeof(sol), :u) && length(sol.u) > 0
    println("First solution element type: ", typeof(sol.u[1]))
    println("First solution element: ", sol.u[1])
end

# Extract trajectory data from solution
# Solution contains 12 elements: [spacecraft_pos(3), spacecraft_mom(3), proxima_pos(3), proxima_mom(3)]
if hasfield(typeof(sol), :u) && length(sol.u) > 0
    t_vals = sol.t
    
    # Extract spacecraft trajectory (first 3 elements are position)
    x_sc = [sol.u[i][1] for i in 1:length(sol.u)]
    y_sc = [sol.u[i][2] for i in 1:length(sol.u)]
    z_sc = [sol.u[i][3] for i in 1:length(sol.u)]
    
    # Extract Proxima trajectory (elements 7-9 are position)
    x_px = [sol.u[i][7] for i in 1:length(sol.u)]
    y_px = [sol.u[i][8] for i in 1:length(sol.u)]
    z_px = [sol.u[i][9] for i in 1:length(sol.u)]
    
    println("Spacecraft trajectory range:")
    println("  X: $(minimum(x_sc)) to $(maximum(x_sc))")
    println("  Y: $(minimum(y_sc)) to $(maximum(y_sc))")
    println("  Z: $(minimum(z_sc)) to $(maximum(z_sc))")
    
    println("Proxima trajectory range:")
    println("  X: $(minimum(x_px)) to $(maximum(x_px))")
    println("  Y: $(minimum(y_px)) to $(maximum(y_px))")
    println("  Z: $(minimum(z_px)) to $(maximum(z_px))")
else
    println("Warning: Could not extract trajectory data from solution")
    x_sc = y_sc = z_sc = [0.0]
    x_px = y_px = z_px = [0.0]
    t_vals = [0.0]
end

# Create 3D trajectory plot for both particles
p3d = plot3d(x_sc, y_sc, z_sc, 
       title="Two-Body Trajectory in Central Scattering",
       xlabel="X Position", 
       ylabel="Y Position", 
       zlabel="Z Position",
       linewidth=2,
       color=:blue,
       label="Spacecraft")

# Add Proxima trajectory
plot3d!(p3d, x_px, y_px, z_px,
        linewidth=2,
        color=:orange,
        label="Proxima Centauri")

# Add initial and final positions for spacecraft
scatter3d!(p3d, [x_sc[1]], [y_sc[1]], [z_sc[1]], 
          color=:green, markersize=6, label="SC Start")
scatter3d!(p3d, [x_sc[end]], [y_sc[end]], [z_sc[end]], 
          color=:red, markersize=6, label="SC End")

# Add initial and final positions for Proxima
scatter3d!(p3d, [x_px[1]], [y_px[1]], [z_px[1]], 
          color=:lightgreen, markersize=6, label="PX Start")
scatter3d!(p3d, [x_px[end]], [y_px[end]], [z_px[end]], 
          color=:darkred, markersize=6, label="PX End")

# Create 2D projection plots for both particles
p1 = plot(x_sc, y_sc, 
          title="XY Projection", 
          xlabel="X Position", 
          ylabel="Y Position",
          linewidth=2, color=:blue, label="Spacecraft")
plot!(p1, x_px, y_px, linewidth=2, color=:orange, label="Proxima")
scatter!(p1, [x_sc[1]], [y_sc[1]], color=:green, markersize=4, label="SC Start")
scatter!(p1, [x_sc[end]], [y_sc[end]], color=:red, markersize=4, label="SC End")
scatter!(p1, [x_px[1]], [y_px[1]], color=:lightgreen, markersize=4, label="PX Start")
scatter!(p1, [x_px[end]], [y_px[end]], color=:darkred, markersize=4, label="PX End")

p2 = plot(x_sc, z_sc, 
          title="XZ Projection", 
          xlabel="X Position", 
          ylabel="Z Position",
          linewidth=2, color=:blue, label="Spacecraft")
plot!(p2, x_px, z_px, linewidth=2, color=:orange, label="Proxima")
scatter!(p2, [x_sc[1]], [z_sc[1]], color=:green, markersize=4, label="SC Start")
scatter!(p2, [x_sc[end]], [z_sc[end]], color=:red, markersize=4, label="SC End")
scatter!(p2, [x_px[1]], [z_px[1]], color=:lightgreen, markersize=4, label="PX Start")
scatter!(p2, [x_px[end]], [z_px[end]], color=:darkred, markersize=4, label="PX End")

p3 = plot(y_sc, z_sc, 
          title="YZ Projection", 
          xlabel="Y Position", 
          ylabel="Z Position",
          linewidth=2, color=:blue, label="Spacecraft")
plot!(p3, y_px, z_px, linewidth=2, color=:orange, label="Proxima")
scatter!(p3, [y_sc[1]], [z_sc[1]], color=:green, markersize=4, label="SC Start")
scatter!(p3, [y_sc[end]], [z_sc[end]], color=:red, markersize=4, label="SC End")
scatter!(p3, [y_px[1]], [z_px[1]], color=:lightgreen, markersize=4, label="PX Start")
scatter!(p3, [y_px[end]], [z_px[end]], color=:darkred, markersize=4, label="PX End")

# Time evolution plot - distances from center of mass
r_sc = sqrt.(x_sc.^2 .+ y_sc.^2 .+ z_sc.^2)
r_px = sqrt.(x_px.^2 .+ y_px.^2 .+ z_px.^2)
r_rel = sqrt.((x_sc .- x_px).^2 .+ (y_sc .- y_px).^2 .+ (z_sc .- z_px).^2)

p4 = plot(t_vals, r_sc,
          title="Distance Evolution vs Time",
          xlabel="Time",
          ylabel="Distance",
          linewidth=2, color=:blue, label="SC from origin")
plot!(p4, t_vals, r_px, linewidth=2, color=:orange, label="PX from origin")
plot!(p4, t_vals, r_rel, linewidth=2, color=:purple, label="Relative distance")

# Display the 3D plot separately
display(p3d)

# Combine 2D projection plots
plot_combined = plot(p1, p2, p3, p4, layout=(2,2), size=(1000,800))

# Display combined plots
display(plot_combined)

println("Trajectory plotting complete!")
println("Spacecraft initial position: ($(x_sc[1]), $(y_sc[1]), $(z_sc[1]))")
println("Spacecraft final position: ($(x_sc[end]), $(y_sc[end]), $(z_sc[end]))")
println("Proxima initial position: ($(x_px[1]), $(y_px[1]), $(z_px[1]))")
println("Proxima final position: ($(x_px[end]), $(y_px[end]), $(z_px[end]))")
println("Minimum relative distance: $(minimum(r_rel))")
println("Maximum relative distance: $(maximum(r_rel))")
