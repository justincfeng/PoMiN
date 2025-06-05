using DoubleFloats

masses = Double64[1, 1.66E-07]        # the Sun and Mercury with mass in units of solar masses

# uncomment for eccentricity = 0
Z_init = Double64[0, 0, 0, 3.8607E+09, 0, 0, 0, -2.67E-12, 0, 0, 2.67E-12, 0]
tspan = (Double64(0), Double64(9.042E+15))
δ = Double64(1.507E+13)

# uncomment for eccentricity = 0.5
# Z_init = Double64[0, 0, 0, 3.8607E+09, 0, 0, 0, -3.27E-12, 0, 0, 3.27E-12, 0]
# tspan = (Double64(0), Double64(1.23E+16))
# δ = Double64(6.15E+13)

# uncomment for eccentricity = 0.9
# Z_init = Double64[0, 0, 0, 3.8607E+09, 0, 0, 0, -3.68E-12, 0, 0, 3.68E-12, 0]
# tspan = (Double64(0), Double64(4.77E+16))
# δ = Double64(4.77E+12)

# uncomment for eccentricity = 0.95
# Z_init = Double64[0, 0, 0, 3.8607E+09, 0, 0, 0, -3.73E-12, 0, 0, 3.73E-12, 0]
# tspan = (Double64(0), Double64(2.70E+17))
# δ = Double64(2.70E+13)

no_adapt = (δ, z, zdot) -> δ                        # turns off adaptive time stepping
adapt = (δ, z, zdot) -> tcour(δ, z, zdot, 0.001)   # last parameter is Courant number
maxit = 30000

sol = hrkintegrator(3, 2, Z_init, x -> pomin.dH(3, masses, x), δ, adapt, tspan, maxit)
#sol = hsintegrator(3, 1, Z_init, x -> pomin.dH(3, masses, x), δ, 1/(δ*100000), tspan, maxit)

println(soln2csv(sol))