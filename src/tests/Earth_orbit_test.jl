using DoubleFloats

masses = Double64[1, 0.00000300336937]        # the Sun and Earth with mass in units of solar masses

Z_init = Double64[0, 0, 0, -1.8667825847E+07, 8.9633743973E+07, 3.8883291466E+07, 0, 0, 0, -2.9839042037E-10, -5.0388889510E-11, -2.1846049069E-11]
tspan = (Double64(0), Double64(1.922E+13))
δ = Double64(1.922E+8)

no_adapt = (δ, z, zdot) -> δ                        # turns off adaptive time stepping
adapt = (δ, z, zdot) -> tcour(δ, z, zdot, 0.0001)   # last parameter is Courant number
maxit = 100000

sol = hrkintegrator(3, 2, Z_init, x -> pomin.dH(3, masses, x), δ, no_adapt, tspan, maxit)
#sol = hsintegrator(3, 1, Z_init, x -> pomin.dH(3, masses, x), δ, 1/(δ*100000), tspan, maxit)

printcsv(sol)