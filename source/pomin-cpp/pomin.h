#pragma once
#include<vector>
#include"particle.h"
#include<string>
#include<math.h>
#include<iostream>
using namespace std;

enum phaseIndex {QX, QY, QZ, PX, PY, PZ};

class Pomin
{
    public:
    
    // Constructor & Destructor
    Pomin(int N);                   // creates a Pomin object with N particles
    ~Pomin();
    
    // functionality
    void read_input(string input_data);

    // computations
    void HamiltonEquations();
    double DHamiltonian(phaseIndex zindex, int c);
    void calcPhi();
    void calcHamiltonian();

    // getters and setters

    int numParts() { return N; };
    void setParticle(int n, double m, double qx, double qy, double qz, double px, double py, double pz);
    int getSmallesta() { return smallesta; };
    int getSmallestb() { return smallestb; };
    double getRsmallest() { return rsmallest; }
    
    private:
    
    
    int N;                          // number of particles
    vector<particle> part;          // N dim array of particles
    double G = 1.0;                 // Newton's Constant (using geometric units)

    int hamiltonian_to_use = 123;   // default is to use all parts of the Hamiltonian

    /* Hamiltonian variables:  H0, H1, H2, H3 from paper
       These are computed by calcHamiltonian() */
    double H[4];

    /* Helper function.  Returns 0, 1.0, or -1.0 */
    double Sign(double x);


    /* PHI variables.  These are computed by calcPhi().
     * These were kept as 1- and 2-dim C-style arrays so that the formula
     * could be easily copied over from the original C code */

    double **r;                     // NxN array of distances between particle pairs
    double *ps;                     // N dim array of momenta squared
    double *E;                      // N dim array of relativistic energies (renamed from mm or En)
    double **y;                     // NxN array of y_ab variables (see paper)
    double **Th;                    // NxN array of Theta_ab variables (see paper)
    double **Xi;                    // NxN array of Xi_ab variables (see paper)

    double *dps;                    // derivative of momenta squared
    double *dE;                     // derivative of relativistic energies
    double **dr;                     // derivative of distances between particle pairs
    double **dy;                    // derivative of y_ab variables (see paper)
    double **dTh;                   // derivative of Theta_ab variabeles (see paper)
    double **dXi;                   // derivative of Xi_ab variables (see paper)

    /* Closest particle pair information.  Computed by calcPhi() */
    double rsmallest;               // smallest distance between any pair of particles
    int smallesta = 0, smallestb = 1;    // particle numbers (a and b) for the closest pair
};
