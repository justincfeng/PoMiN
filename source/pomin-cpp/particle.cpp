#include "particle.h"

particle::particle(double mass, double qx, double qy, double qz, double px, double py, double pz)
{
    m = mass;
    q[0] = qx;
    q[1] = qy;
    q[2] = qz;
    p[0] = px;
    p[1] = py;
    p[2] = pz;
}