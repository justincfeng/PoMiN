############################################################
# eccentricity.py - calculate orbital eccentricity
############################################################
# This file will calculate the orbital momentum required for a given
# eccentricity. The values must be hard-coded below.
############################################################
#!/usr/bin/python

# sympy allows for symbolic manipulations.
from sympy import *

# Define a variable, p, for momentum.
p = Symbol('p')
T = Symbol('T')

# Set the values for a desired orbit.
G=1                    # gravitational constant
M=1                    # mass of central object
m=1.66e-7              # mass of orbiting object
r=3.8607e9             # radius at perihthesis
e=0.95                 # eccentricity of the orbit

# Define the equations to use.
r_ap = r*((1+e)/(1-e)) # radius at apthesis
b = sqrt(r*r_ap)       # semi-minor axis
a = b/(sqrt(1-e**2))   # semi-major axis

k = G*M*m              # spring constant
L = r*p                # angular momentum
E = (p**2/(2*m)) - k/r # energy

# Solve for the orbital period.
print(solve(a - ((G*M*T**2)/(4*pi**2))**(1/3), T))
# Solve the eccentricity equation for momentum.
print(solve(e - sqrt(1 + (2*E*L**2)/(k**2*m)), p))

############################################################
# END eccentricity.py
############################################################

