using DoubleFloats

masses = Double64[1, 0.00000300336937, 0.1221, 1.0788, 0.9092]        # Sun, Earth, Proxima, Alpha A, and Alpha B masses in units of solar masses

Z_init = Double64[0, 0, 0, -1.8667825847E+07, 8.9633743973E+07, 3.8883291466E+07, -9.9157696012E+12, -7.5891180437E+12, -2.4171267294E+13, -1.0514396567E+13, -8.7920419041E+12, -2.4558076843E+13, -1.0514148949E+13, -8.7899759136E+12, -2.4558922297E+13, 0, 0, 0, -2.9839042037E-10, -5.0388889510E-11, -2.1846049069E-11, -3.8172314445E-06, 9.0494364315E-06, 8.9932012183E-06, -3.2681576591E-05, 8.2784394612E-05, 7.2543810790E-05, -3.5191057704E-05, 7.8741673392E-05, 7.9956568487E-05]
tspan = (Double64(0), Double64(1.359815e+14)) # 21.22 years
δ = Double64(9.61E+9) # about half a day

no_adapt = (δ, z, zdot) -> δ                        # turns off adaptive time stepping
adapt = (δ, z, zdot) -> tcour(δ, z, zdot, 0.0001)   # last parameter is Courant number
maxit = 100000

sol = hrkintegrator(3, 5, Z_init, x -> pomin.dH(3, masses, x), δ, no_adapt, tspan, maxit)
#sol = hsintegrator(3, 1, Z_init, x -> pomin.dH(3, masses, x), δ, 1/(δ*100000), tspan, maxit)

printcsv(sol)