using DoubleFloats

masses = Double64[1]        # one object, the Sun, with mass in units of solar masses

Z_init = Double64[1.7552537847E+17 0 0 0 8.0722511038E-04 0]
tspan = (Double64(0), Double64(1.5E+21))
δ = Double64(1.5E+17)
tadap = t -> 1.0  # effectively turns off adaptive time stepping since δ will be multiplied by 1 each timestep
maxit = 100000

sol = hrkintegrator(3, 1, Z_init, x -> pomin.dH(3, masses, x), δ, tadap, tspan, maxit)
#sol = hsintegrator(3, 1, Z_init, x -> pomin.dH(3, masses, x), δ, 1/(δ*100000), tspan, maxit)

println(soln2csv(sol))