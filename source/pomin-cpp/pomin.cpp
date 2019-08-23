#include "pomin.h"

Pomin::Pomin(int n) : N(n)
{
    // memory allocation

    r = new double*[N];
    ps = new double[N];
    E = new double[N];
    y = new double*[N];
    Th = new double*[N];
    Xi = new double*[N];
    for (int i=0; i<N; i++)
    {
        r[i] = new double[N];
        y[i] = new double[N];
        Th[i] = new double[N];
        Xi[i] = new double[N];
    }
}

Pomin::~Pomin()
{
}

void Pomin::setParticle(int n, double m, double qx, double qy, double qz, double px, double py, double pz)
{
    if (n < N) {
        part[n].mass = m;
        part[n].qx = qx;
    }
}


void Pomin::read_input(string input_data)
{
    istringstream input(input_data);
    
    string description;
    double start_time, end_time, timestep, grav, courant = 1.0;
    int interations;
    

    
    // Populate pomin object with input data
    double m, qx, qy, qz, px, py, pz;
    for(int i=0; i<N; i++)
    {
        input >> m >> "," >> qx >> "," >>
    }
    
}

void Pomin::HamiltonEquations()
{
    /************************************************************
     * COMPUTE HAMILTON EQUATIONS
     ************************************************************/
    for(int n=0; n<N; n++) {
        for (int x=0; x<3; x++) {
            part[n].pdot[x] = -1.0 * DHamiltonian(static_cast<phaseIndex>(x),n);
            part[n].qdot[x] = DHamiltonian(static_cast<phaseIndex>(x+3),n);
        }
    }

}

double Pomin::DHamiltonian(phaseIndex zindex, int c)
{
    /************************************************************
     * DHamiltonian
     ************************************************************/
    /* DHamiltonian calculates the derivative of the Hamiltonian with respect
     * to a given phase space coordinate. */
    /************************************************************
     * zindex is the phase space index (0=qx,1=qy,2=qz,3=px,4=py,5=pz).
     * c      is the particle index.
     *
     * qk & pk are arrays where k is a spatial index (0=x,1=y,2=z).
     * dps, dE, dr, dy, dTh, dXi are derivatives of ps, mm, m, r, y, Th, Xi. 
     ************************************************************/


    /************************************************************
     * CALCULATE DERIVATIVES OF PHI
     ************************************************************
     * The simplest way to do this is to write a nested loop over the
     * particle number. However, this is not the most efficient way to
     * do it. This function (DHamiltonian) is called in a loop within
     * the function "HamiltonEquations," which runs from 1 to N, where N
     * is the number of particles. If we have nested loops over particle
     * labels within in this function, the number of computations will
     * scale as N^3.
     *
     * Here, we take advantage of the fact that the only nonvanishing
     * derivatives come from terms with particle labels that match the
     * particle label for the variable that we're taking the derivative
     * with respect to (in our case, the label c).
     *
     * When implementing this, we set one of the particle labels to c,
     * and sum over the remaining particle label. We do separate sums
     * for both a and b to ensure that all the nonvanishing terms are
     * calculated.
     *
     * This allows us to avoid nested loops in this function--since this
     * function contains loops over particle number, the number of
     * computations scale as N^2.
     ************************************************************/
    if(zindex<3) {
		
		// For Test particles
        if(c>=N) {
		dps[c] = 0;
        dE[c] = 0;
		}
		
        for(int a=0; a<N; a++) {
            dps[a] = 0;
            dE[a] = 0;

            if(a!=c) {
                dr[a][c]  = (1.0/r[a][c]) * (part[c].q[zindex]-part[a].q[zindex]);
                dTh[a][c] = (1.0/r[a][c]) * (part[a].p[zindex]-Th[a][c]*dr[a][c]);
                dXi[a][c] = 0;

                if(part[a].m > 0) {
                    dy[a][c] = (E[a]*dTh[a][c]*Th[a][c] + dE[a]*(-part[a].m*part[a].m - Th[a][c]*Th[a][c]))/(E[a]*E[a]*sqrt(part[a].m*part[a].m + Th[a][c]*Th[a][c]));
                }
                else {
                    dy[a][c] = Sign(Th[a][c])*(E[a]*dTh[a][c] - dE[a]*Th[a][c])/(E[a]*E[a]);
                }
            }
        }
        for(int b=0; b<N; b++) {
            if(b!=c) {
                dr[c][b]  = (1.0/r[c][b]) * (part[b].q[zindex]-part[c].q[zindex]) * (-1.0);
                dTh[c][b] = (1.0/r[c][b]) * (part[c].p[zindex]*(-1.0) - Th[c][b]*dr[c][b]);
                dXi[c][b] = 0;

                if(part[c].m > 0) {
                    dy[c][b] = (E[c]*dTh[c][b]*Th[c][b] + dE[c]*(-part[c].m*part[c].m - Th[c][b]*Th[c][b]))/(E[c]*E[c]*sqrt(part[c].m*part[c].m + Th[c][b]*Th[c][b]));
                }
                else {
                    dy[c][b] = Sign(Th[c][b])*(E[c]*dTh[c][b] - dE[c]*Th[c][b])/(E[c]*E[c]);
                }
            }
        }
    }

    else if(zindex > 2) {
        // This allows us to use zindex directly as a fixed subscript.
        zindex = static_cast<phaseIndex>(zindex - 3);

        for(int a=0; a<N; a++) {
            if(a==c) {
                dps[a] = 2.0*part[a].p[zindex];
            }
            else {
                dps[a] = 0;
            }

            dE[a] = (0.5/(E[a])) * dps[a];
        }
        
            // For Test particles
            if(c>=N) {
			dps[c] = 2.0*part[c].p[zindex];
            dE[c] = (0.5/(E[c])) * dps[c];
			}

        for(int a=0; a<N; a++) {
            if(a!=c) {
                dr[a][c]  = 0;
                dTh[a][c] = 0;
                dXi[a][c] = part[a].p[zindex];

                if(part[a].m > 0) {
                    dy[a][c] = (E[a]*dTh[a][c]*Th[a][c] + dE[a]*(-part[a].m*part[a].m - Th[a][c]*Th[a][c]))/(E[a]*E[a]*sqrt(part[a].m*part[a].m + Th[a][c]*Th[a][c]));
                }
                else {
                    dy[a][c] = Sign(Th[a][c])*(E[a]*dTh[a][c] - dE[a]*Th[a][c])/(E[a]*E[a]);
                }
            }
        }

        for(int b=0; b<N; b++) {
            if(b!=c) {
                dr[c][b]  = 0;
                dTh[c][b] = (1.0/r[c][b]) * (part[b].q[zindex]-part[c].q[zindex]);
                dXi[c][b] = part[b].p[zindex];

                if(part[c].m > 0) {
                    dy[c][b] = (E[c]*dTh[c][b]*Th[c][b] + dE[c]*(-part[c].m*part[c].m - Th[c][b]*Th[c][b]))/(E[c]*E[c]*sqrt(part[c].m*part[c].m + Th[c][b]*Th[c][b]));
                }
                else {
                    dy[c][b] = Sign(Th[c][b])*(E[c]*dTh[c][b] - dE[c]*Th[c][b])/(E[c]*E[c]);
                }
            }
        }
    }

    /************************************************************
     * CALCULATE dH (the derivative of the hamiltonian)
     ************************************************************/
    // dH[4] = {dHnewt,dH1,dH2,dH3}
    double dH[4] = {0,0,0,0};

    for(int n=0; n<N; n++)  {
        dH[1] += dE[n];
        dH[0] += 0.5*dps[n]/part[n].m;
        
        // For Test particles
        if(c>=N) {
		dH[1] += dE[c];
        dH[0] += 0.5*dps[c]/part[c].m;
		}
    }

    for(int a=0; a<N; a++) {
        if(c!=a) {
            dH[1] -= (G*0.5)*((-(E[a]*E[c]*(E[c]*E[c]*ps[a] + E[a]*E[a]*(E[c]*E[c] + ps[c]))*dr[a][c]) + (-(dE[a]*E[c]*E[c]*E[c]*ps[a]) + E[a]*E[c]*E[c]*(dps[a]*E[c] + dE[c]*ps[a]) + E[a]*E[a]*E[a]*(dps[c]*E[c] + dE[c]*(E[c]*E[c] - ps[c])) + dE[a]*E[a]*E[a]*E[c]*(E[c]*E[c] + ps[c]))*r[a][c])/(E[a]*E[a]*E[c]*E[c]*r[a][c]*r[a][c]));
            dH[2] += (G*0.25)*((7.0*dXi[a][c]*r[a][c] - dTh[c][a]*r[a][c]*Th[a][c] - dTh[a][c]*r[a][c]*Th[c][a] + dr[a][c]*Th[a][c]*Th[c][a] - 7.0*dr[a][c]*Xi[a][c])/(r[a][c]*r[a][c]));

            if(part[c].m > 0) {
                dH[3] -= (G*0.25)*(-((2.0*E[a]*E[c]*dy[c][a]*r[a][c]*y[c][a]*(2.0*ps[c]*ps[c]*Th[a][c]*Th[a][c] + 4.0*ps[c]*Th[a][c]*Th[c][a]*Xi[a][c] - 2.0*(ps[c] - 2.0*Th[c][a]*Th[c][a])*Xi[a][c]*Xi[a][c] + E[c]*E[c]*(2.0*(-(Th[a][c]*Th[c][a]) + Xi[a][c])*(-(Th[a][c]*Th[c][a]) + Xi[a][c]) + Th[a][c]*Th[c][a]*(Th[a][c]*Th[c][a] - 8.0*Xi[a][c])*y[c][a] - ps[a]*Th[c][a]*Th[c][a]*(2.0 + 3.0*y[c][a]) + ps[c]*(ps[a]*y[c][a] - Th[a][c]*Th[a][c]*(2.0 + 3.0*y[c][a])))) + dE[c]*E[a]*r[a][c]*y[c][a]*(1.0 + y[c][a])*(2.0*ps[c]*ps[c]*Th[a][c]*Th[a][c] + 4.0*ps[c]*Th[a][c]*Th[c][a]*Xi[a][c] - 2.0*(ps[c] - 2.0*Th[c][a]*Th[c][a])*Xi[a][c]*Xi[a][c] + E[c]*E[c]*(2.0*(-(Th[a][c]*Th[c][a]) + Xi[a][c])*(-(Th[a][c]*Th[c][a]) + Xi[a][c]) + Th[a][c]*Th[c][a]*(Th[a][c]*Th[c][a] - 8.0*Xi[a][c])*y[c][a] - ps[a]*Th[c][a]*Th[c][a]*(2.0 + 3.0*y[c][a]) + ps[c]*(ps[a]*y[c][a] - Th[a][c]*Th[a][c]*(2.0 + 3.0*y[c][a])))) - E[a]*E[c]*dy[c][a]*r[a][c]*(1.0 + y[c][a])*(-2.0*ps[c]*ps[c]*Th[a][c]*Th[a][c] - 4.0*Th[c][a]*Th[c][a]*Xi[a][c]*Xi[a][c] + 2.0*ps[c]*Xi[a][c]*(-2.0*Th[a][c]*Th[c][a] + Xi[a][c]) + E[c]*E[c]*(-2.0*(-(Th[a][c]*Th[c][a]) + Xi[a][c])*(-(Th[a][c]*Th[c][a]) + Xi[a][c]) + Th[a][c]*Th[c][a]*(-(Th[a][c]*Th[c][a]) + 8.0*Xi[a][c])*y[c][a] + ps[a]*Th[c][a]*Th[c][a]*(2.0 + 3.0*y[c][a]) + ps[c]*(-(ps[a]*y[c][a]) + Th[a][c]*Th[a][c]*(2.0 + 3.0*y[c][a])))) - E[a]*E[c]*dr[a][c]*y[c][a]*(1.0 + y[c][a])*(-2.0*ps[c]*ps[c]*Th[a][c]*Th[a][c] - 4.0*Th[c][a]*Th[c][a]*Xi[a][c]*Xi[a][c] + 2.0*ps[c]*Xi[a][c]*(-2.0*Th[a][c]*Th[c][a] + Xi[a][c]) + E[c]*E[c]*(-2.0*(-(Th[a][c]*Th[c][a]) + Xi[a][c])*(-(Th[a][c]*Th[c][a]) + Xi[a][c]) + Th[a][c]*Th[c][a]*(-(Th[a][c]*Th[c][a]) + 8.0*Xi[a][c])*y[c][a] + ps[a]*Th[c][a]*Th[c][a]*(2.0 + 3.0*y[c][a]) + ps[c]*(-(ps[a]*y[c][a]) + Th[a][c]*Th[a][c]*(2.0 + 3.0*y[c][a])))) - dE[a]*E[c]*r[a][c]*y[c][a]*(1.0 + y[c][a])*(-2.0*ps[c]*ps[c]*Th[a][c]*Th[a][c] - 4.0*Th[c][a]*Th[c][a]*Xi[a][c]*Xi[a][c] + 2.0*ps[c]*Xi[a][c]*(-2.0*Th[a][c]*Th[c][a] + Xi[a][c]) + E[c]*E[c]*(-2.0*(-(Th[a][c]*Th[c][a]) + Xi[a][c])*(-(Th[a][c]*Th[c][a]) + Xi[a][c]) + Th[a][c]*Th[c][a]*(-(Th[a][c]*Th[c][a]) + 8.0*Xi[a][c])*y[c][a] + ps[a]*Th[c][a]*Th[c][a]*(2.0 + 3.0*y[c][a]) + ps[c]*(-(ps[a]*y[c][a]) + Th[a][c]*Th[a][c]*(2.0 + 3.0*y[c][a])))) + E[a]*r[a][c]*y[c][a]*(1.0 + y[c][a])*(2.0*dps[c]*E[c]*E[c]*E[c]*Th[a][c]*Th[a][c] + 4.0*ps[c]*ps[c]*Th[a][c]*(-(E[c]*dTh[a][c]) + dE[c]*Th[a][c]) - 4.0*E[c]*E[c]*E[c]*dXi[a][c]*Xi[a][c] + 4.0*E[c]*E[c]*E[c]*dTh[c][a]*Th[a][c]*Xi[a][c] + 2.0*dps[c]*E[c]*Xi[a][c]*Xi[a][c] - dps[c]*E[c]*E[c]*E[c]*ps[a]*y[c][a] + 3.0*dps[c]*E[c]*E[c]*E[c]*Th[a][c]*Th[a][c]*y[c][a] + 8.0*E[c]*E[c]*E[c]*dTh[c][a]*Th[a][c]*Xi[a][c]*y[c][a] + Th[c][a]*Th[c][a]*(-8.0*E[c]*dXi[a][c]*Xi[a][c] + 8.0*dE[c]*Xi[a][c]*Xi[a][c] + E[c]*E[c]*E[c]*(dy[c][a]*(3.0*ps[a] - Th[a][c]*Th[a][c]) - 2.0*dTh[a][c]*Th[a][c]*(2.0 + y[c][a]) + dps[a]*(2.0 + 3.0*y[c][a]))) + ps[c]*((-4.0*dps[c]*E[c] + 3.0*E[c]*E[c]*E[c]*dy[c][a])*Th[a][c]*Th[a][c] + 4.0*E[c]*(dXi[a][c] - dTh[a][c]*Th[c][a])*Xi[a][c] - 4.0*dE[c]*Xi[a][c]*Xi[a][c] - E[c]*E[c]*E[c]*(ps[a]*dy[c][a] + dps[a]*y[c][a]) + Th[a][c]*(-4.0*E[c]*dTh[c][a]*Xi[a][c] + Th[c][a]*(-4.0*E[c]*dXi[a][c] + 8.0*dE[c]*Xi[a][c]) + 2.0*E[c]*E[c]*E[c]*dTh[a][c]*(2.0 + 3.0*y[c][a]))) + Th[c][a]*(-4.0*E[c]*Xi[a][c]*(dps[c]*Th[a][c] + 2.0*dTh[c][a]*Xi[a][c]) + E[c]*E[c]*E[c]*(4.0*dXi[a][c]*Th[a][c]*(1.0 + 2.0*y[c][a]) + 4.0*Xi[a][c]*(2.0*dy[c][a]*Th[a][c] + dTh[a][c]*(1.0 + 2.0*y[c][a])) + 2.0*dTh[c][a]*(-(Th[a][c]*Th[a][c]*(2.0 + y[c][a])) + ps[a]*(2.0 + 3.0*y[c][a]))))))/(E[a]*E[a]*E[c]*E[c]*E[c]*E[c]*r[a][c]*r[a][c]*y[c][a]*y[c][a]*(1.0 + y[c][a])*(1.0 + y[c][a])*(1.0 + y[c][a]))));
                dH[0] += 0.5*G*part[a].m*part[c].m*dr[a][c]/(r[a][c]*r[a][c]);
            }
            else {
                double signTh = Sign(Th[c][a]);
                dH[3] -= (G*0.25)*(2.0*E[a]*ps[c]*dr[a][c]*(1.0 + y[c][a])*(8.0*signTh*sqrt(ps[c])*Th[a][c]*Xi[a][c]*y[c][a] - 4.0*Xi[a][c]*Xi[a][c]*y[c][a] - ps[c]*Th[a][c]*Th[a][c]*(-1.0 + y[c][a])*(3.0 + y[c][a]) +  ps[a]*ps[c]*(1.0 + y[c][a])*(-1.0 + 3.0*y[c][a])) + r[a][c]*(dps[c]*E[a]*(1.0 + y[c][a])*(8.0*signTh*sqrt(ps[c])*Th[a][c]*Xi[a][c]*y[c][a] - 2.0*6.0*Xi[a][c]*Xi[a][c]*y[c][a] - ps[c]*Th[a][c]*Th[a][c]*(3.0 + y[c][a]*(2.0 + y[c][a])) + ps[a]*ps[c]*(1.0 + y[c][a]*(2.0 + 3.0*y[c][a]))) - 2.0*(dps[a]*E[a]*ps[c]*ps[c]*(1.0 + y[c][a])*(1.0 + y[c][a])*(-1.0 + 3.0*y[c][a]) + dE[a]*ps[c]*(1.0 + y[c][a])*(-8.0*signTh*sqrt(ps[c])*Th[a][c]*Xi[a][c]*y[c][a] + 4.0*Xi[a][c]*Xi[a][c]*y[c][a] + ps[c]*Th[a][c]*Th[a][c]*(-1.0 + y[c][a])*(3.0 + y[c][a]) - ps[a]*ps[c]*(1.0 + y[c][a])*(-1.0 + 3.0*y[c][a])) + 2.0*E[a]*sqrt(ps[c])*(dTh[c][a]*(1.0 + y[c][a])*(4.0*sqrt(ps[c])*Th[a][c]*Xi[a][c] - 4.0*signTh*Xi[a][c]*Xi[a][c] - signTh*ps[c]*Th[a][c]*Th[a][c]*(2.0 + y[c][a]) + signTh*ps[a]*ps[c]*(2.0 + 3.0*y[c][a])) + sqrt(ps[c])*(4.0*dXi[a][c]*(signTh*sqrt(ps[c])*Th[a][c] - Xi[a][c])*y[c][a]*(1.0 + y[c][a]) - sqrt(ps[c])*dTh[a][c]*(1.0 + y[c][a])*(-4.0*signTh*Xi[a][c]*y[c][a] + sqrt(ps[c])*Th[a][c]*(-1.0 + y[c][a])*(3.0 + y[c][a])) + dy[c][a]*(-8.0*signTh*sqrt(ps[c])*Th[a][c]*Xi[a][c]*y[c][a] + 2.0*Xi[a][c]*Xi[a][c]*(1.0 + 3.0*y[c][a]) + ps[c]*(-3.0*ps[a]*y[c][a]*(1.0 + y[c][a]) + Th[a][c]*Th[a][c]*(-2.0 + y[c][a]*(3.0 + y[c][a])))))))))/(2.0*E[a]*E[a]*sqrt(ps[c])*ps[c]*r[a][c]*r[a][c]*(1.0 + y[c][a])*(1.0 + y[c][a])*(1.0 + y[c][a]));
            }
        }
    }
    for(int b=0; b<N; b++) {
        if(b!=c) {
            dH[1] -= (G*0.5)*((-(E[c]*E[b]*(E[b]*E[b]*ps[c] + E[c]*E[c]*(E[b]*E[b] + ps[b]))*dr[c][b]) + (-(dE[c]*E[b]*E[b]*E[b]*ps[c]) + E[c]*E[b]*E[b]*(dps[c]*E[b] + dE[b]*ps[c]) + E[c]*E[c]*E[c]*(dps[b]*E[b] + dE[b]*(E[b]*E[b] - ps[b])) + dE[c]*E[c]*E[c]*E[b]*(E[b]*E[b] + ps[b]))*r[c][b])/(E[c]*E[c]*E[b]*E[b]*r[c][b]*r[c][b]));
            dH[2] += (G*0.25)*((7.0*dXi[c][b]*r[c][b] - dTh[b][c]*r[c][b]*Th[c][b] - dTh[c][b]*r[c][b]*Th[b][c] + dr[c][b]*Th[c][b]*Th[b][c] - 7.0*dr[c][b]*Xi[c][b])/(r[c][b]*r[c][b]));

            if(part[b].m > 0) {
                dH[3] -= (G*0.25)*(-((2.0*E[c]*E[b]*dy[b][c]*r[c][b]*y[b][c]*(2.0*ps[b]*ps[b]*Th[c][b]*Th[c][b] + 4.0*ps[b]*Th[c][b]*Th[b][c]*Xi[c][b] - 2.0*(ps[b] - 2.0*Th[b][c]*Th[b][c])*Xi[c][b]*Xi[c][b] + E[b]*E[b]*(2.0*(-(Th[c][b]*Th[b][c]) + Xi[c][b])*(-(Th[c][b]*Th[b][c]) + Xi[c][b]) + Th[c][b]*Th[b][c]*(Th[c][b]*Th[b][c] - 8.0*Xi[c][b])*y[b][c] - ps[c]*Th[b][c]*Th[b][c]*(2.0 + 3.0*y[b][c]) + ps[b]*(ps[c]*y[b][c] - Th[c][b]*Th[c][b]*(2.0 + 3.0*y[b][c])))) + dE[b]*E[c]*r[c][b]*y[b][c]*(1.0 + y[b][c])*(2.0*ps[b]*ps[b]*Th[c][b]*Th[c][b] + 4.0*ps[b]*Th[c][b]*Th[b][c]*Xi[c][b] - 2.0*(ps[b] - 2.0*Th[b][c]*Th[b][c])*Xi[c][b]*Xi[c][b] + E[b]*E[b]*(2.0*(-(Th[c][b]*Th[b][c]) + Xi[c][b])*(-(Th[c][b]*Th[b][c]) + Xi[c][b]) + Th[c][b]*Th[b][c]*(Th[c][b]*Th[b][c] - 8.0*Xi[c][b])*y[b][c] - ps[c]*Th[b][c]*Th[b][c]*(2.0 + 3.0*y[b][c]) + ps[b]*(ps[c]*y[b][c] - Th[c][b]*Th[c][b]*(2.0 + 3.0*y[b][c])))) - E[c]*E[b]*dy[b][c]*r[c][b]*(1.0 + y[b][c])*(-2.0*ps[b]*ps[b]*Th[c][b]*Th[c][b] - 4.0*Th[b][c]*Th[b][c]*Xi[c][b]*Xi[c][b] + 2.0*ps[b]*Xi[c][b]*(-2.0*Th[c][b]*Th[b][c] + Xi[c][b]) + E[b]*E[b]*(-2.0*(-(Th[c][b]*Th[b][c]) + Xi[c][b])*(-(Th[c][b]*Th[b][c]) + Xi[c][b]) + Th[c][b]*Th[b][c]*(-(Th[c][b]*Th[b][c]) + 8.0*Xi[c][b])*y[b][c] + ps[c]*Th[b][c]*Th[b][c]*(2.0 + 3.0*y[b][c]) + ps[b]*(-(ps[c]*y[b][c]) + Th[c][b]*Th[c][b]*(2.0 + 3.0*y[b][c])))) - E[c]*E[b]*dr[c][b]*y[b][c]*(1.0 + y[b][c])*(-2.0*ps[b]*ps[b]*Th[c][b]*Th[c][b] - 4.0*Th[b][c]*Th[b][c]*Xi[c][b]*Xi[c][b] + 2.0*ps[b]*Xi[c][b]*(-2.0*Th[c][b]*Th[b][c] + Xi[c][b]) + E[b]*E[b]*(-2.0*(-(Th[c][b]*Th[b][c]) + Xi[c][b])*(-(Th[c][b]*Th[b][c]) + Xi[c][b]) + Th[c][b]*Th[b][c]*(-(Th[c][b]*Th[b][c]) + 8.0*Xi[c][b])*y[b][c] + ps[c]*Th[b][c]*Th[b][c]*(2.0 + 3.0*y[b][c]) + ps[b]*(-(ps[c]*y[b][c]) + Th[c][b]*Th[c][b]*(2.0 + 3.0*y[b][c])))) - dE[c]*E[b]*r[c][b]*y[b][c]*(1.0 + y[b][c])*(-2.0*ps[b]*ps[b]*Th[c][b]*Th[c][b] - 4.0*Th[b][c]*Th[b][c]*Xi[c][b]*Xi[c][b] + 2.0*ps[b]*Xi[c][b]*(-2.0*Th[c][b]*Th[b][c] + Xi[c][b]) + E[b]*E[b]*(-2.0*(-(Th[c][b]*Th[b][c]) + Xi[c][b])*(-(Th[c][b]*Th[b][c]) + Xi[c][b]) + Th[c][b]*Th[b][c]*(-(Th[c][b]*Th[b][c]) + 8.0*Xi[c][b])*y[b][c] + ps[c]*Th[b][c]*Th[b][c]*(2.0 + 3.0*y[b][c]) + ps[b]*(-(ps[c]*y[b][c]) + Th[c][b]*Th[c][b]*(2.0 + 3.0*y[b][c])))) + E[c]*r[c][b]*y[b][c]*(1.0 + y[b][c])*(2.0*dps[b]*E[b]*E[b]*E[b]*Th[c][b]*Th[c][b] + 4.0*ps[b]*ps[b]*Th[c][b]*(-(E[b]*dTh[c][b]) + dE[b]*Th[c][b]) - 4.0*E[b]*E[b]*E[b]*dXi[c][b]*Xi[c][b] + 4.0*E[b]*E[b]*E[b]*dTh[b][c]*Th[c][b]*Xi[c][b] + 2.0*dps[b]*E[b]*Xi[c][b]*Xi[c][b] - dps[b]*E[b]*E[b]*E[b]*ps[c]*y[b][c] + 3.0*dps[b]*E[b]*E[b]*E[b]*Th[c][b]*Th[c][b]*y[b][c] + 8.0*E[b]*E[b]*E[b]*dTh[b][c]*Th[c][b]*Xi[c][b]*y[b][c] + Th[b][c]*Th[b][c]*(-8.0*E[b]*dXi[c][b]*Xi[c][b] + 8.0*dE[b]*Xi[c][b]*Xi[c][b] + E[b]*E[b]*E[b]*(dy[b][c]*(3.0*ps[c] - Th[c][b]*Th[c][b]) - 2.0*dTh[c][b]*Th[c][b]*(2.0 + y[b][c]) + dps[c]*(2.0 + 3.0*y[b][c]))) + ps[b]*((-4.0*dps[b]*E[b] + 3.0*E[b]*E[b]*E[b]*dy[b][c])*Th[c][b]*Th[c][b] + 4.0*E[b]*(dXi[c][b] - dTh[c][b]*Th[b][c])*Xi[c][b] - 4.0*dE[b]*Xi[c][b]*Xi[c][b] - E[b]*E[b]*E[b]*(ps[c]*dy[b][c] + dps[c]*y[b][c]) + Th[c][b]*(-4.0*E[b]*dTh[b][c]*Xi[c][b] + Th[b][c]*(-4.0*E[b]*dXi[c][b] + 8.0*dE[b]*Xi[c][b]) + 2.0*E[b]*E[b]*E[b]*dTh[c][b]*(2.0 + 3.0*y[b][c]))) + Th[b][c]*(-4.0*E[b]*Xi[c][b]*(dps[b]*Th[c][b] + 2.0*dTh[b][c]*Xi[c][b]) + E[b]*E[b]*E[b]*(4.0*dXi[c][b]*Th[c][b]*(1.0 + 2.0*y[b][c]) + 4.0*Xi[c][b]*(2.0*dy[b][c]*Th[c][b] + dTh[c][b]*(1.0 + 2.0*y[b][c])) + 2.0*dTh[b][c]*(-(Th[c][b]*Th[c][b]*(2.0 + y[b][c])) + ps[c]*(2.0 + 3.0*y[b][c]))))))/(E[c]*E[c]*E[b]*E[b]*E[b]*E[b]*r[c][b]*r[c][b]*y[b][c]*y[b][c]*(1.0 + y[b][c])*(1.0 + y[b][c])*(1.0 + y[b][c]))));
                dH[0] += 0.5*G*part[c].m*part[b].m*dr[c][b]/(r[c][b]*r[c][b]);
            }
            else {
                double signTh = Sign(Th[b][c]);
                dH[3] -= (G*0.25)*(2.0*E[c]*ps[b]*dr[c][b]*(1.0 + y[b][c])*(8.0*signTh*sqrt(ps[b])*Th[c][b]*Xi[c][b]*y[b][c] - 4.0*Xi[c][b]*Xi[c][b]*y[b][c] - ps[b]*Th[c][b]*Th[c][b]*(-1.0 + y[b][c])*(3.0 + y[b][c]) +  ps[c]*ps[b]*(1.0 + y[b][c])*(-1.0 + 3.0*y[b][c])) + r[c][b]*(dps[b]*E[c]*(1.0 + y[b][c])*(8.0*signTh*sqrt(ps[b])*Th[c][b]*Xi[c][b]*y[b][c] - 2.0*6.0*Xi[c][b]*Xi[c][b]*y[b][c] - ps[b]*Th[c][b]*Th[c][b]*(3.0 + y[b][c]*(2.0 + y[b][c])) + ps[c]*ps[b]*(1.0 + y[b][c]*(2.0 + 3.0*y[b][c]))) - 2.0*(dps[c]*E[c]*ps[b]*ps[b]*(1.0 + y[b][c])*(1.0 + y[b][c])*(-1.0 + 3.0*y[b][c]) + dE[c]*ps[b]*(1.0 + y[b][c])*(-8.0*signTh*sqrt(ps[b])*Th[c][b]*Xi[c][b]*y[b][c] + 4.0*Xi[c][b]*Xi[c][b]*y[b][c] + ps[b]*Th[c][b]*Th[c][b]*(-1.0 + y[b][c])*(3.0 + y[b][c]) - ps[c]*ps[b]*(1.0 + y[b][c])*(-1.0 + 3.0*y[b][c])) + 2.0*E[c]*sqrt(ps[b])*(dTh[b][c]*(1.0 + y[b][c])*(4.0*sqrt(ps[b])*Th[c][b]*Xi[c][b] - 4.0*signTh*Xi[c][b]*Xi[c][b] - signTh*ps[b]*Th[c][b]*Th[c][b]*(2.0 + y[b][c]) + signTh*ps[c]*ps[b]*(2.0 + 3.0*y[b][c])) + sqrt(ps[b])*(4.0*dXi[c][b]*(signTh*sqrt(ps[b])*Th[c][b] - Xi[c][b])*y[b][c]*(1.0 + y[b][c]) - sqrt(ps[b])*dTh[c][b]*(1.0 + y[b][c])*(-4.0*signTh*Xi[c][b]*y[b][c] + sqrt(ps[b])*Th[c][b]*(-1.0 + y[b][c])*(3.0 + y[b][c])) + dy[b][c]*(-8.0*signTh*sqrt(ps[b])*Th[c][b]*Xi[c][b]*y[b][c] + 2.0*Xi[c][b]*Xi[c][b]*(1.0 + 3.0*y[b][c]) + ps[b]*(-3.0*ps[c]*y[b][c]*(1.0 + y[b][c]) + Th[c][b]*Th[c][b]*(-2.0 + y[b][c]*(3.0 + y[b][c])))))))))/(2.0*E[c]*E[c]*sqrt(ps[b])*ps[b]*r[c][b]*r[c][b]*(1.0 + y[b][c])*(1.0 + y[b][c])*(1.0 + y[b][c]));
            }
        }
    }


    switch (hamiltonian_to_use) {
        case 0:         dH[1]=dH[2]=dH[3]=0; break;
        case 1:   dH[0]      =dH[2]=dH[3]=0; break;
        case 12:  dH[0]            =dH[3]=0; break;
        case 13:  dH[0]      =dH[2]      =0; break;
        default: case 123: dH[0]                  =0; break;
    }

    return(dH[0]+dH[1]+dH[2]+dH[3]);
}

void Pomin::calcPhi()
// precondition: Pomin has been populated with q and p data
// postcondition: all of the "Phi" variables r, ps, E, y, Th, Xi have been calculated
{
    /************************************************************
     * COMPUTE: ps, E
     ************************************************************/
    for(int n=0; n<N; n++) {
        ps[n] = part[n].p[0]*part[n].p[0] + part[n].p[1]*part[n].p[1] + part[n].p[2]*part[n].p[2];
        E[n] = sqrt(part[n].m*part[n].m + ps[n]);

        if(E[n] == 0) {
            cerr << "Pomin: 0 energy for particle %i.\n" << n+1 << endl;
            exit(EXIT_FAILURE);
        }
    }

    /************************************************************
     * COMPUTE: r, Th, y, Xi
     ************************************************************/
    // Start with the assumption that the first two particles are closest.
    // ats & bts are the labels for the closest particles.

    for (int a=0; a<N; a++) {
        for (int b=0; b<N; b++) {
            if (b != a) {
                r[a][b] = sqrt((part[b].q[0]-part[a].q[0])*(part[b].q[0]-part[a].q[0]) 
                        + (part[b].q[1]-part[a].q[1])*(part[b].q[1]-part[a].q[1]) 
                        + (part[b].q[2]-part[a].q[2])*(part[b].q[2]-part[a].q[2]));

                if(r[a][b] == 0) {
                    cerr << "Pomin: 0 separation between particles " << a+1 << " & " << b+1 << endl;
                    exit(EXIT_FAILURE);
                }

                // This sets the value of rsmallest in the first iteration of the loop.
                rsmallest=r[smallesta][smallestb];

                // Check if the current pair of particles are closer.
                if(r[a][b] < rsmallest) {
                    rsmallest = r[a][b];
                    smallesta = a;
                    smallestb = b;
                }

                Th[b][a] = (-1.0/r[a][b]) * (part[b].p[0] * (part[b].q[0]-part[a].q[0]) \
                        + part[b].p[1] * (part[b].q[1]-part[a].q[1]) \
                        + part[b].p[2] * (part[b].q[2]-part[a].q[2]));
                y[b][a]  = (1.0/E[b]) * sqrt(part[b].m*part[b].m + Th[b][a]*Th[b][a]);
                Xi[a][b] = part[a].p[0]*part[b].p[0] \
                           + part[a].p[1]*part[b].p[1] \
                           + part[a].p[2]*part[b].p[2];
            }
        }
    }
}

void Pomin::calcHamiltonian()
// precondition: calcPhi() has been called
// postcondition: The Hamiltonian variables H0, H1, H2, H3 have been calculated
{
   

    /************************************************************
     * HAMILTONIAN CALCULATOR
     ************************************************************/
    // Does not include contributions from test particles.
    // Subtract the rest mass from H1 to magnify changes in H.
    for(int n=0; n<N; n++) {
        H[1] += E[n] - part[n].m;

        if(part[n].m > 0) {
            H[0] += 0.5*ps[n]/part[n].m;
        }
    }
/*
    for(int a=0; a<N; a++) {
        for(int b=0; b<N; b++) {
            if(b != a) {
                H[1] -= (G*0.5)*(E[a]*E[b]/r[a][b])*(1.0 + (ps[a]/(E[a]*E[a])) + (ps[b]/(E[b]*E[b])));
                H[2] += (G*.25)*(((-1.0)*(Th[a][b]*Th[b][a]) + 7.0*Xi[a][b])/r[a][b]);

                if(part[b].m > 0) {
                    H[3] -= (G*.25)*(1.0/r[a][b])*(1.0/(E[a]*E[b]*(y[b][a]+1.0)*(y[b][a]+1.0)*y[b][a]))*( 2.0*( 2.0*Xi[a][b]*Xi[a][b]*Th[b][a]*Th[b][a] - 2.0*(Th[a][b])*(-Th[b][a])*(Xi[a][b])*ps[b] + Th[a][b]*Th[a][b]*ps[b]*ps[b] - Xi[a][b]*Xi[a][b]*ps[b] )*(1.0/(E[b]*E[b])) + 2.0*(-ps[a]*Th[b][a]*Th[b][a] + Th[a][b]*Th[a][b]*Th[b][a]*Th[b][a] + 2.0*(Th[a][b])*(-Th[b][a])*(Xi[a][b]) + Xi[a][b]*Xi[a][b] - Th[a][b]*Th[a][b]*ps[b]) + ((-3.0*ps[a]*Th[b][a]*Th[b][a] + Th[a][b]*Th[a][b]*Th[b][a]*Th[b][a] + 8.0*(Th[a][b])*(-Th[b][a])*(Xi[a][b]) + ps[a]*ps[b] - 3.0*Th[a][b]*Th[a][b]*ps[b] )*y[b][a] ) );
                    H[0] -= 0.5*G*part[a].m*part[b].m/r[a][b];
                }
                else {
                    H[3] -= (G*0.25)*(ps[b]*Th[a][b]*Th[a][b]*(-1.0 + y[b][a])*(3.0 + y[b][a]) - ps[a]*ps[b]*(1.0 + y[b][a])*(-1.0 + 3.0*y[b][a]) + 4.0*Xi[a][b]*(-2.0*Th[a][b]*Th[b][a] + Xi[a][b]*y[b][a]))/(sqrt(ps[a])*sqrt(ps[b])*r[a][b]*(1.0 + y[b][a])*(1.0 + y[b][a]));
                }
            }
        }
    }

    switch (hamiltonian_to_use) {
        case 0:        H[1]=H[2]=H[3]=0; break;
        case 1:   H[0]     =H[2]=H[3]=0; break;
        case 12:  H[0]          =H[3]=0; break;
        case 13:  H[0]     =H[2]     =0; break;
        default: case 123: H[0]               =0; break;
    }
*/
}

double Pomin::Sign(double x) {
    /* Return the sign of A. */
    return (x > 0) ? 1.0 : ((x < 0) ? -1.0 : 0);
}
