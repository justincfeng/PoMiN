using DoubleFloats

# Mass of Sun, Earth, Moon 
masses = Double64[1, 3.0033693739E-06, 3.6950832738E-08]

# Initial Z: Positions first, momenta second
Z_init = [0, 0, 0, 99615330, 0, 0, 99888800, 0, 24626.7, 0, 0, 0, 0, 3.03440999999999E-10, 0, 0, 3.85106199999999E-12, 0]


tspan = (Double64(0), Double64(4.9997E+12)) 
# tspan = (Double64(0), Double64(1.922E+11)) # 10 days
δ = Double64(1.7E+8)

no_adapt = (δ, z, zdot) -> δ                        # turns off adaptive time stepping
adapt = (δ, z, zdot) -> tcour(δ, z, zdot, 0.0001)   # last parameter is Courant number
maxit = 100000

sol = hrkintegrator(3, 3, Z_init, x -> pomin.dH(3, masses, x), δ, no_adapt, tspan, maxit)

printcsv(sol)
