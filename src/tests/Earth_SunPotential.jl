using DoubleFloats

# Mass of Earth 
masses = Double64[3.0033693739E-06]
# masses = Double64[1, 3.0033693739E-06]

# Initial Z: Positions first, momenta second
Z_init = [99615330, 0, 0, 0, 3.03440999999999E-10, 0]
# Z_init = [0,0,0, 99615330, 0, 0, 0,0,0, 0, 3.03440999999999E-10, 0]

tspan = (Double64(0), Double64(6.4E+13))  # about 10 years
δ = Double64(1.75E+9)

no_adapt = (δ, z, zdot) -> δ                        # turns off adaptive time stepping
adapt = (δ, z, zdot) -> tcour(δ, z, zdot, 0.0001)   # last parameter is Courant number
maxit = 100000

sol = hrkintegrator(3, 1, Z_init, x -> pomin.dH_plus_Sun(3, masses, x), δ, no_adapt, tspan, maxit)
# sol = hrkintegrator(3, 2, Z_init, x -> pomin.dH(3, masses, x), δ, no_adapt, tspan, maxit)

printcsv(sol)
