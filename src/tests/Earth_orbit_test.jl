using DoubleFloats

masses = Double64[1, 0.00000300336937]        # the Sun and Earth with mass in units of solar masses

Z_init = Double64[0, 0, 0, -1.866782584728260E+07, 8.96337439732320E+07, 3.888329146643650E+07, 0, 0, 0, -2.983904203718850E-10, -5.038888951024780E-11, -2.184604906902130E-11]  # ICRS coords J2000, init conditions from astropy
tspan = (Double64(0), Double64(6.407E+13)) # about 10 years
δ = Double64(6.407E+8) # 1e5 time steps

no_adapt = (δ, z, zdot) -> δ                        # turns off adaptive time stepping
adapt = (δ, z, zdot) -> tcour(δ, z, zdot, 0.0001)   # last parameter is Courant number
maxit = 100000

sol = hrkintegrator(3, 2, Z_init, x -> pomin.dH(3, masses, x), δ, no_adapt, tspan, maxit)
#sol = hsintegrator(3, 1, Z_init, x -> pomin.dH(3, masses, x), δ, 1/(δ*100000), tspan, maxit)

printcsv(sol)