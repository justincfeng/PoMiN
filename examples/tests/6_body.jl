using DoubleFloats

masses = Double64[1, 3.0033693739E-06, 1.0057832537E-33, 0.1221, 1.0788, 0.9092] # Sun, Earth, Spacechip, Proxima, Alpha A, Alpha B in M

Z_init = Double64[0, 0, 0, -1.8667826140E+07, 8.9633743993E+07, 3.8883291714E+07, -1.8667826140E+07, 8.9894055560E+07, 3.8883291714E+07, -9.9157694997E+12, -7.5891180191E+12, -2.4171267258E+13, -1.0514396564E+13, -8.7920419314E+12, -2.4558076740E+13, -1.0514148951E+13, -8.7899759755E+12, -2.4558922285E+13, 0, 0, 0, -2.9839042080E-10, -5.0388889700E-11, -2.1846049145E-11, -7.1894181893E-35, -5.4924334992E-35, -1.7509992075E-34, -3.8172314441E-06, 9.0494364315E-06, 8.9932012189E-06, -3.2681576581E-05, 8.2784394610E-05, 7.2543810801E-05, -2.9658611095E-05, 6.6362559735E-05, 6.7386459079E-05]
# Position first: Sun, Earth, Spacechip, Proxima, Alpha A, Alpha B
# Momentum second: Sun, Earth, Spacechip, Proxima, Alpha A, Alpha B


tspan = (Double64(0), Double64(1.922E+14)) # 30 years
#tspan = (Double64(0), Double64(1.36039051E+14)) # 21.2325 years
δ = Double64(1.922E+10) # about one day

no_adapt = (δ, z, zdot) -> δ                        # turns off adaptive time stepping
adapt = (δ, z, zdot) -> tcour(δ, z, zdot, 0.0001)   # last parameter is Courant number
maxit = 100000

sol = hrkintegrator(3, 5, Z_init, x -> pomin.dH(3, masses, x), δ, no_adapt, tspan, maxit)
#sol = hsintegrator(3, 1, Z_init, x -> pomin.dH(3, masses, x), δ, 1/(δ*100000), tspan, maxit)

printcsv(sol)