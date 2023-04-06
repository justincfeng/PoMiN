using DoubleFloats

masses = Double64[1]        # one object, the Sun, with mass in units of solar masses

Z_init = Double64[1.7552537847E+17 0 0 0 8.0722511038E-04 0]
tspan = (Double64(0), Double64(7.5E+21))
δ = Double64(1.5E+17)
noadapt = (δ,z,zdot) -> δ  # effectively turns off adaptive time stepping
maxit = 100000

sol = hrkintegrator(3, 1, Z_init, x -> pomin.dH_plus_MW(3, masses, x), δ, noadapt, tspan, maxit)
#sol = hsintegrator(3, 1, Z_init, x -> pomin.dH(3, masses, x), δ, 1/(δ*100000), tspan, maxit)

printcsv(sol)