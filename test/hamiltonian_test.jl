include("../source/HamPM.jl")

n = 4;
d = 3;

m = ones(n)      ;
z = zeros(n*2*d) ;

i=3;

z[2*d*(i-1)+1] = 1. ; 
z[2*d*(i-1)+2] = 1. ; 
z[2*d*(i-1)+3] = 1. ; 

z[2*d*(i-1)+1+n] = 2. ; 
z[2*d*(i-1)+2+n] = 2. ; 
z[2*d*(i-1)+3+n] = 2. ; 

@test Z2q( n , d , i , z ) == ones(d)

@test Z2p( n , d , i , z ) == 2*ones(d)
