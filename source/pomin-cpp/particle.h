#pragma once

class particle
{
    public:
        // constructor
        particle(double m, double qx, double qy, double qz, double px, double py, double pz);

        // position and its time derivative
        double q[3];
        double qdot[3];
        // momentum and its time derivative
        double p[3];
        double pdot[3];
        // mass
        double m;

        bool isTestPart;
};