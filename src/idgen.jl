

m1 = 1.0 ;
m2 = 0.000000166 ;

q1 = [  0. ;  0. ;  0. ] ;
q2 = [  1. ;  0. ;  0. ]*3860700000 ;

pmag = 2.6716223e-12 ;

p1 = [  0. ; -1. ;  0. ]*pmag ;
p2 = [  0. ;  1. ;  0. ]*pmag ;

part2b = pomin.Particles([m1,m2],[q1,q2],[p1,p2]) ;

mvec = pomin.Part2m(part2b) ;
Zvec = pomin.Part2Z(part2b) ;

ω = 0.01 ;

δ = 1.507e13 ;

tspan = (0.,1.507e15) ;

dHf = z-> pomin.dH( 3 , mvec , z )

pomin.hsintegrator( Zvec , dHf , δ , ω , tspan , 100 )

# CHECK DIMENSIONS

