using DoubleFloats

masses = Double64[1, 0.00000300336937, 0.1221]        # Sun, Earth, and Proxima masses in units of solar masses

Z_init = Double64[0, 0, 0, -1.8667825847E+07, 8.9633743973E+07, 3.8883291466E+07, -9.9157696012E+12, -7.5891180437E+12, -2.4171267294E+13, 0, 0, 0, -2.9839042037E-10, -5.0388889510E-11, -2.1846049069E-11, -3.8172314445E-06, 9.0494364315E-06, 8.9932012183E-06]
tspan = (Double64(0), Double64(1.922E+14))
δ = Double64(1.922E+10)

no_adapt = (δ, z, zdot) -> δ                        # turns off adaptive time stepping
adapt = (δ, z, zdot) -> tcour(δ, z, zdot, 0.0001)   # last parameter is Courant number
maxit = 10000

sol = hrkintegrator(3, 3, Z_init, x -> pomin.dH(3, masses, x), δ, no_adapt, tspan, maxit)
#sol = hsintegrator(3, 1, Z_init, x -> pomin.dH(3, masses, x), δ, 1/(δ*100000), tspan, maxit)

printcsv(sol)