/************************************************************ 
 * PoMiN: A Hamiltonian Post-Minkowski N-Body Code
 *
 * See pomin(1) for copyright, authorship, and usage info.
 ************************************************************/

/************************************************************
 * HEADERS
 ************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <unistd.h>
#include <time.h>
#include <string.h>

// sys/time has millisecond resolution.
#include <sys/time.h>
// This macro in required to print the uint65_t type.
#define __STDC_FORMAT_MACROS
// timespec struct needs uint65_t, provided by inttypes.
#include <inttypes.h>

/************************************************************
 * SET MATHEMATICAL PRECISION
 ************************************************************/
/* There is no floating-point literal type-specifier for the quadmath
 * library; i.e., we cannot write "1.0q" and expect this 1.0 to be of
 * the __float128 type.  Therefore, we use preprocessor definitions to
 * force literals of this type so that [quad x 1.0 = quad] and not
 * [quad x 1.0 = double] since double is the default type for
 * floating-point literals. */

#ifdef QUAD
#include <quadmath.h>
#define PRECISION      __float128
#define pow            powq
#define sqrt           sqrtq
#define signbit        signbitq
#define NEGATIVEONE    strtoflt128("-1.0", NULL)
#define ZERO           strtoflt128("0.0", NULL)
#define QUARTER        strtoflt128("0.25", NULL)
#define HALF           strtoflt128("0.5", NULL)
#define ONE            strtoflt128("1.0", NULL)
#define TWO            strtoflt128("2.0", NULL)
#define THREE          strtoflt128("3.0", NULL)
#define FOUR           strtoflt128("4.0", NULL)
#define FIVE           strtoflt128("5.0", NULL)
#define SIX            strtoflt128("6.0", NULL)
#define SEVEN          strtoflt128("7.0", NULL)
#define EIGHT          strtoflt128("8.0", NULL)
#define NINE           strtoflt128("9.0", NULL)
#define qsixth         ONE/SIX
#define SIXTH          qsixth

#else // not QUAD
#include <math.h>
#define PRECISION      double
#define NEGATIVEONE    -1.0
#define ZERO           0.0
#define QUARTER        0.25
#define HALF           0.5
#define ONE            1.0
#define TWO            2.0
#define THREE          3.0
#define FOUR           4.0
#define FIVE           5.0
#define SIX            6.0
#define SEVEN          7.0
#define EIGHT          8.0
#define NINE           9.0
#define dsixth         ONE/SIX
#define SIXTH          dsixth
#endif // QUAD

/************************************************************
 * GLOBAL VARIABLES
 ************************************************************/
// N is the number of particles.
static unsigned int N, Ntot, hamiltonian_to_use = 123;
static PRECISION gravity, timestep, maxtimestep, *mass;
static char *name;

// The following 2D arrays are indexed by: [particle #][dimension].
// q    & p    are the canonical coordinates (q) and conjugate momenta (p).
// qdot & pdot are their time derivatives.
// Q & P are intermediate values in the RK4 algorithm.
PRECISION (*q)[3], (*qdot)[3], (*Q)[3];
PRECISION (*p)[3], (*pdot)[3], (*P)[3];

PRECISION *ps, *En, **r, **y, **Th, **Xi;
PRECISION *dps, *dEn, **dr, **dy, **dTh, **dXi;

/************************************************************
 * FUNCTION PROTOTYPES
 ************************************************************/
// Computes the right hand side of Hamilton's equations.
void HamiltonEquations(PRECISION *H, PRECISION courant, bool calculate_hamiltonian);

// Computes the derivatives of the Hamiltonian.
PRECISION DHamiltonian(unsigned int zindex, unsigned const int c);

/* Change the boolean output of signbit to the __float128/double valued
 * outputs, 1.0 or -1.0 */
PRECISION Sign(PRECISION A);

/************************************************************
 * MAIN
 ************************************************************/
// Contains Runge-Kutta integration loop.
int main(int argc, char **argv) {

    /************************************************************
     * PARSE COMMAND LINE ARGUMENTS
     ************************************************************/
    // Set default values.
    bool error = false, debug = false, verbose = false;
    // Save the name of the program into a variable to identify output.
    name = argv[0];

    // Loop through the coEnand line arguments.
    for (int i=1; i<argc; i++)  {

        // Check if argument is flagged.
        if (!strncmp(argv[i], "-", 1)) {

            if (!strcmp(argv[i], "-d")) {
                debug = true;
            }
            else if (!strcmp(argv[i], "-v")) {
                verbose = true;
            }
            else if (!strcmp(argv[i], "-p")) {
                hamiltonian_to_use = strtol(argv[++i], NULL, 10);

                if (hamiltonian_to_use != 0  && hamiltonian_to_use != 1  \
                        && hamiltonian_to_use != 12 && hamiltonian_to_use != 13 \
                        && hamiltonian_to_use != 123)
                {
                    error = true;
                }
            }
            else { // unrecognized flag
		fprintf(stderr, "unrecognized flag %s\n", argv[i]);
                error = true;
            }
        } 
        else { // non-flagged argument
	    fprintf(stderr, "non-flagged argument %s\n", argv[i]);
            error = true;
        }
    }

    // Fail if run without input.
    if (isatty(fileno(stdin))) {
        error = true;
    }

    if (error) {
        fprintf(stderr, "usage: %s [-d] [-p <part>] \
                \ncalculate motions of an n-body system\n\n \
                input\t\t\tis taken from stdin (see pomin(1) for more)\n \
                output\t\t\tis given to stdout\n \
                -d\t\t\tprint debug data to stdout\n \
                -p {1,12,13,123}\tuse parts 1, 12, 13, or 123 of the Hamiltonian\n \
                \ne.g., grep <unique-id> <input-file> | %s > <output-file>\n", name, name);
        exit(EXIT_FAILURE);
    }

    /************************************************************
     * SAVE GIT COMMIT HASH & STATE
     ************************************************************/
    /* To aid in the project development, it would be useful to know
     * which version of the code is associated with any particular set
     * of output values; this will save the git hash and state (dirty
     * vs. clean) of the working directory. */
    char githash[43], gitstate[6];

    // Open processes to obtain said information.
    FILE *hash  = popen("git rev-parse --verify HEAD", "r");
    FILE *state = popen("git status --porcelain",      "r");

    // Save output of those coEnands to variables.
    fgets(githash,  sizeof(githash)-1,  hash);
    fgets(gitstate, sizeof(gitstate)-1, state);

    // Close the processes.
    pclose(hash);
    pclose(state);

    /* Set state of working directory (i.e., if gitstate is empty, set
     * it to "clean"). */
    snprintf(gitstate, sizeof(gitstate), ((!strlen(gitstate)) ? "clean" : "dirty"));
    // The ternary conditional a?b:c executes b if a is true or c otherwise.

    /************************************************************
     * GET CURRENT TIME AND DATE
     ************************************************************/
    char datetime[20];
    struct tm *sTm;
    time_t now = time(0);

    sTm = gmtime(&now);
    // Formatted in accordance with ISO 8601 (YYYY-MM-DDTHH:MM:SS).
    strftime(datetime, sizeof(datetime), "%Y-%m-%dT%H:%M:%S", sTm);

    /************************************************************
     * READ LINE FROM STDIN
     ************************************************************/
    char *line = NULL;
    size_t size;

    // TODO: use recursion with a different (not main) function.
    if (getline(&line, &size, stdin) == -1) {
        exit(EXIT_SUCCESS);
    }

    // Determine the number of particles given.
    for (int i=0; line[i]!='\0'; i++) {
        // i.e., count the number of coEnas in the line.
        if (line[i] == ',') {
            N++;

            // End if there's nothing between two coEnas (an empty entry).
            if(line[i+1] == ',') {
                break;
            }
        }
    }
    // Account for preceding (7) & per-particle (7) entries.
    N-=7; N/=7;
    
    Ntot = N;

    // Disregard the description (used only to identify input).
    char* description = strtok(line, ",");

    /************************************************************
     * MEMORY ALLOCATION (main)
     ************************************************************/
    PRECISION *H    = malloc(sizeof(PRECISION[4]));
    mass            = malloc(sizeof(PRECISION[N]));
    q               = malloc(sizeof(PRECISION[N][3]));
    p               = malloc(sizeof(PRECISION[N][3]));
    qdot            = malloc(sizeof(PRECISION[N][3]));
    pdot            = malloc(sizeof(PRECISION[N][3]));

    /* Q & P are the phase space coordinates used
     * to compute the right hand side of Hamilton's Equations at 
     * each partial timestep in the RK4 scheme. */
    Q = malloc(sizeof(PRECISION[N][3]));
    P = malloc(sizeof(PRECISION[N][3]));
    PRECISION (*qi)[3] = malloc(sizeof(PRECISION[N][3]));
    PRECISION (*pi)[3] = malloc(sizeof(PRECISION[N][3]));

    // Exit if memory has been depleted.
    if (H == NULL || mass == NULL ||  q == NULL || p == NULL || \
            qdot == NULL || pdot == NULL || Q == NULL || P == NULL || \
            pi == NULL || qi == NULL)
    {
        fprintf(stderr, "%s: out of memory\n", name);
        exit(EXIT_FAILURE);
    }

    ps  	= malloc(sizeof(PRECISION[N]));
    En  	= malloc(sizeof(PRECISION[N]));
    r  		= malloc(N * sizeof(PRECISION *));
    y  		= malloc(N * sizeof(PRECISION *));
    Th 		= malloc(N * sizeof(PRECISION *));
    Xi 		= malloc(N * sizeof(PRECISION *));
    dps  	= malloc(sizeof(PRECISION[N]));
    dEn  	= malloc(sizeof(PRECISION[N]));
    dr  	= malloc(N * sizeof(PRECISION *));
    dy  	= malloc(N * sizeof(PRECISION *));
    dTh 	= malloc(N * sizeof(PRECISION *));
    dXi 	= malloc(N * sizeof(PRECISION *));

    /* Quit if out of memory. */
    if(ps == NULL || En == NULL || r  == NULL ||      \
            y  == NULL || Th == NULL || Xi == NULL)
    {
        fprintf(stderr, "%s: out of memory\n", name);
        exit(EXIT_FAILURE);
    }
    if(dps == NULL || dEn == NULL || \
            dr == NULL || dy == NULL || dTh == NULL || dXi == NULL)
    {
        fprintf(stderr, "%s: out of memory\n", name);
        exit(EXIT_FAILURE);
    }

    for(int n=0; n<N; n++) {
        r[n]  = malloc(N * sizeof(PRECISION));
        y[n]  = malloc(N * sizeof(PRECISION));
        Th[n] = malloc(N * sizeof(PRECISION));
        Xi[n] = malloc(N * sizeof(PRECISION));

        if(r[n] == NULL || y[n] == NULL || Th[n] == NULL || Xi[n] == NULL) {
            fprintf(stderr, "%s: out of memory\n", name);
            exit(EXIT_FAILURE);
        }
    }

    for(int n=0; n<N; n++)  {
        dr[n]  = malloc(N * sizeof(PRECISION));
        dy[n]  = malloc(N * sizeof(PRECISION));
        dTh[n] = malloc(N * sizeof(PRECISION));
        dXi[n] = malloc(N * sizeof(PRECISION));

        if(dr[n] == NULL || dy[n] == NULL || dTh[n] == NULL || dXi[n] == NULL)
        {
            fprintf(stderr, "%s: out of memory\n", name);
            exit(EXIT_FAILURE);
        }
    }

#ifdef QUAD
    /************************************************************
     * CREATE QUAD PRECISION VALUES
     ************************************************************/
    /* In order to use the quadmath.h library, numbers must be read in
     * as strings, then converted to __float128 type with the
     * strtoflt128 function, and finally converted back to strings with
     * quadmath_snprintf for printing. */

    // Create strings to read and print quad precision numbers.
    char    time_str[512], duration_str[512],     timestep_str[512];
    char courant_str[512],  gravity_str[512];
    char mass_str[N][512];

    char    H_str[3][512];
    char    TotH_str[512];
    char    q_str[N][3][512],    p_str[N][3][512];
    char    qdot_str[N][3][512], pdot_str[N][3][512];

    // Read some values from input file.
    snprintf(time_str,        sizeof(time_str),        "%s", strtok(NULL, ","));
    snprintf(duration_str,    sizeof(duration_str),    "%s", strtok(NULL, ","));
    snprintf(timestep_str,    sizeof(timestep_str),    "%s", strtok(NULL, ","));
    unsigned const int iterations         = strtol(strtok(NULL, ","), NULL, 10);
    snprintf(courant_str,     sizeof(courant_str),     "%s", strtok(NULL, ","));
    snprintf(gravity_str,     sizeof(gravity_str),     "%s", strtok(NULL, ","));
    unsigned const int numspecparticles   = strtol(strtok(NULL, ","), NULL, 10);

    // Convert strings to quad numbers.
    PRECISION time        = strtoflt128(time_str,        NULL);
    PRECISION duration    = strtoflt128(duration_str,    NULL);
    timestep              = strtoflt128(timestep_str,    NULL);
    PRECISION courant     = strtoflt128(courant_str,     NULL);
    gravity               = strtoflt128(gravity_str,     NULL);

    // Save quad numbers to strings for printing.
    quadmath_snprintf(time_str,        sizeof(time_str),        "%.35Qf", time);
    quadmath_snprintf(duration_str,    sizeof(duration_str),    "%.35Qf", duration);
    quadmath_snprintf(timestep_str,    sizeof(timestep_str),    "%.35Qf", timestep);
    quadmath_snprintf(courant_str,     sizeof(courant_str),     "%.35Qf", courant);
    quadmath_snprintf(gravity_str,     sizeof(gravity_str),     "%.35Qf", gravity);

    // Now do all of the above for each of the particles.
    for(int n=0; n<N; n++) {
        snprintf(mass_str[n], sizeof(mass_str[n]), "%s", strtok(NULL, ","));
        mass[n] = strtoflt128(mass_str[n],   NULL);
        quadmath_snprintf(mass_str[n], sizeof(mass_str[n]),  "%.35Qf", mass[n]);

        for(int x=0; x<3; x++) {
            snprintf(q_str[n][x], sizeof(q_str[n][x]), "%s", strtok(NULL, ","));
            q[n][x] = strtoflt128(q_str[n][x], NULL);
            quadmath_snprintf(q_str[n][x], sizeof(q_str[n][x]),"%.35Qf", q[n][x]);
        }

        for(int x=0; x<3; x++) {
            snprintf(p_str[n][x], sizeof(p_str[n][x]), "%s", strtok(NULL, ","));
            p[n][x] = strtoflt128(p_str[n][x], NULL);
            quadmath_snprintf(p_str[n][x], sizeof(p_str[n][x]),"%.35Qf", p[n][x]);
        }
    }

#else /* not QUAD */

    PRECISION time        = strtod(strtok(NULL, ","), NULL);
    PRECISION duration    = strtod(strtok(NULL, ","), NULL);
    timestep              = strtod(strtok(NULL, ","), NULL);
    unsigned const int iterations       = strtol(strtok(NULL, ","), NULL, 10);
    PRECISION courant     = strtod(strtok(NULL, ","), NULL);
    gravity               = strtod(strtok(NULL, ","), NULL);
    unsigned const int numspecparticles = strtol(strtok(NULL, ","), NULL, 10);

    for(int n=0; n<N; n++) {
        mass[n] = strtod(strtok(NULL, ","), NULL);
        if(mass[n] == 0 && hamiltonian_to_use == 0){
            fprintf(stderr, "%s: zero mass with Newtonian Hamiltonian\n", name);
            exit(EXIT_FAILURE);
        }
        for (int x=0; x<3; x++) {
            q[n][x] = strtod(strtok(NULL, ","), NULL);
        }
        for (int x=0; x<3; x++) {
            p[n][x] = strtod(strtok(NULL, ","), NULL);
        }
    }

#endif // QUAD

    // Initialize qi and pi to q and p.
    memcpy(qi,q,sizeof(PRECISION[N][3]));
    memcpy(pi,p,sizeof(PRECISION[N][3]));

    /************************************************************
     * PRINT METADATA & INPUT
     ************************************************************/
    if (debug) {
#if QUAD
        // Grab precision to print as output below.
        char* precision = "QUAD";
#else // not QUAD
        char* precision = "DOUBLE";
#endif // QUAD

        printf("YYYY-MM-DDTHH:MM:SS (ISO 8601) of %s code using H:%i with Git hash:\n"\
                "%s [%s] %s\ndescription,start-time,end-time,timestep,iterations,"\
                "courant-number,gravitational-constant",\
                precision, hamiltonian_to_use, datetime, gitstate, githash);

        for(int n=1; n<N+1; n++) {
            printf(",mass_%i,qx_%i,qy_%i,qz_%i,px_%i,py_%i,pz_%i", n,n,n,n,n,n,n);
        }

        printf("\n%s", description);

#ifdef QUAD
        printf(",%s,%s,%s,%i,%s,%s", time_str, duration_str, timestep_str, \
                iterations, courant_str, gravity_str);

        for(int n=0; n<N; n++) {
            printf(",%s,%s,%s,%s,%s,%s,%s", mass_str[n], q_str[n][0], q_str[n][1], q_str[n][2],\
                    p_str[n][0], p_str[n][1], p_str[n][2]);
        }

#else // not QUAD
        printf(",%.35lf,%lf,%lf,%i,%lf,%lf", time, duration, timestep, \
                iterations, courant, gravity);

        for(int n=0; n<N; n++) {
            printf(",%.35lf,%.35lf,%.35lf,%.35lf,%.35lf,%.35lf,%.35lf", \
                    mass[n], q[n][0], q[n][1], q[n][2], p[n][0], p[n][1], p[n][2]);
        }
#endif // QUAD
        printf("\n\n");
    }

    // Free the input buffer.
    free(line);

    /************************************************************
     * PRINT HEADER
     ************************************************************/
    printf("iteration,time");

    for(int n=1; n<N+1; n++) {

        printf(",qx_%i,qy_%i,qz_%i,px_%i,py_%i,pz_%i", n,n,n,n,n,n);
    
	if (debug) printf(",qdotx_%i,qdoty_%i,qdotz_%i,pdotx_%i,pdoty_%i,pdotz_%i", n,n,n,n,n,n);
    }

    if (debug) printf(",H0,H1,H2,H3,Hamiltonian");
	
    printf("\n");

    /************************************************************
     * PRINT INITIAL DATA
     ************************************************************/
#ifdef QUAD
    printf("0,%s", time_str);

    for(int n=0; n<N; n++) {
        printf(",%s,%s,%s,%s,%s,%s", q_str[n][0], q_str[n][1], q_str[n][2],\
                p_str[n][0], p_str[n][1], p_str[n][2]);
    }

#else // not QUAD
    printf("0,%.35lf", time);

    for(int n=0; n<N; n++) {
        for(int x=0; x<3; x++) {
            printf(",%.35lf", q[n][x]);
        }
        for(int x=0; x<3; x++) {
            printf(",%.35lf", p[n][x]);
        }
    }
    
#endif // QUAD

	if(numspecparticles>0) {
    N = numspecparticles;
	}

/*
    if (debug) {
	printf(",0,0,0,0,0,0,0");
    }
  */      
    printf("\n");

    // Wall time in milliseconds (wall_time = end-start).
    struct timespec start, end;

    if (verbose) {
        // Record the start time.
#ifndef DARWIN
        clock_gettime(CLOCK_MONOTONIC_RAW, &start);
#endif
    }

    // Set Maximum Timestep
    maxtimestep=timestep;

    /************************************************************
     * MAIN LOOP
     ************************************************************/
    for(int iteration=1; time<duration; iteration++) {
        /************************************************************
         * Zero Out Arrays
         ************************************************************/
        //memset(ps,0,N*sizeof(PRECISION));
        //memset(En,0,N*sizeof(PRECISION));
        //memset(r,0,N*sizeof(*PRECISION));		
        //memset(y,0,N*N*sizeof(PRECISION));
        //memset(Th,0,N*N*sizeof(PRECISION));
        //memset(Xi,0,N*N*sizeof(PRECISION));

        //memset(dps,0,N*sizeof(PRECISION));
        //memset(dEn,0,N*sizeof(PRECISION));
        //memset(dr,0,N*N*sizeof(PRECISION));
        //memset(dy,0,N*N*sizeof(PRECISION));
        //memset(dTh,0,N*N*sizeof(PRECISION));
        //memset(dXi,0,N*N*sizeof(PRECISION));


        /************************************************************
         * RK4 FIRST STAGE
         ************************************************************/
        /* Each stage of the fourth order Runge-Kutta integration scheme
         * involves the calculation of the quantities: kp,kq.
         *
         * These quantities are the right-hand side of Hamilton's
         * Equations, computed with the function HamiltonEquations.
         ************************************************************/

        // Initialize Q and P to qi and pi.
        memcpy(Q,qi,sizeof(PRECISION[Ntot][3]));
        memcpy(P,pi,sizeof(PRECISION[Ntot][3]));

        HamiltonEquations(H, courant, true);
        time += timestep;

        /************************************************************
         * RK4 SECOND STAGE
         ************************************************************/
        for(int n=0; n<Ntot; n++) {
            for(int x=0; x<3; x++) {
                /* This prepares the next set of phase space variables for
                 * re-calculating the quantities kp,kq needed for the second stage
                 * of the Runge-Kutta algorithm. */
                Q[n][x] = qi[n][x]+HALF*timestep*qdot[n][x];
                P[n][x] = pi[n][x]+HALF*timestep*pdot[n][x];

                /* The updated values for the p's and q's are obtained by computing
                 * a weighted average of pdot and qdot. The following calculation
                 * produces the first piece of the weighted average. */
                q[n][x] = qi[n][x]+(SIXTH)*timestep*qdot[n][x];
                p[n][x] = pi[n][x]+(SIXTH)*timestep*pdot[n][x];
            }
        }

        HamiltonEquations(H, 0, false);

        /************************************************************
         * RK4 THIRD STAGE: COMPUTING pdot3,qdot3
         ************************************************************/
        for(int n=0; n<Ntot; n++) {
            for(int x=0; x<3; x++) {
                Q[n][x] = qi[n][x]+HALF*timestep*qdot[n][x];
                P[n][x] = pi[n][x]+HALF*timestep*pdot[n][x];

                q[n][x] += (SIXTH)*timestep*(TWO*qdot[n][x]);
                p[n][x] += (SIXTH)*timestep*(TWO*pdot[n][x]);
            }
        }

        HamiltonEquations(H, 0, false);

        /************************************************************
         * RK4 FOURTH STAGE: COMPUTING pdot4,qdot4
         ************************************************************/
        for(int n=0; n<Ntot; n++) {
            for(int x=0; x<3; x++) {
                Q[n][x] = qi[n][x]+timestep*qdot[n][x];
                P[n][x] = pi[n][x]+timestep*pdot[n][x];

                q[n][x] += (SIXTH)*timestep*(TWO*qdot[n][x]);
                p[n][x] += (SIXTH)*timestep*(TWO*pdot[n][x]);
            }
        }

        HamiltonEquations(H, 0, false);

        /************************************************************
         * RK4 FINISHING & PRINTING OUTPUT
         ************************************************************/
        printf("%i", iteration);
#ifdef QUAD
        quadmath_snprintf(time_str, 512, "%.35Qf", time);
        printf(",%s", time_str);
#else // not QUAD
        printf(",%.35lf", time);
#endif

        for(int n=0; n<Ntot; n++) {
            for(int x=0; x<3; x++) {
                qi[n][x] = q[n][x] += (SIXTH)*timestep*qdot[n][x];
                pi[n][x] = p[n][x] += (SIXTH)*timestep*pdot[n][x];
            }

#ifdef QUAD
            for(int x=0; x<3; x++) {
                quadmath_snprintf(q_str[n][x], 512, "%.35Qf", q[n][x]);
                printf(",%s", q_str[n][x]);
            }
            for(int x=0; x<3; x++) {
                quadmath_snprintf(p_str[n][x], 512, "%.35Qf", p[n][x]);
                printf(",%s", p_str[n][x]);
            }

            if (debug) {

                for(int x=0; x<3; x++) {
                    quadmath_snprintf(qdot_str[n][x], 512, "%.35Qf", qdot[n][x]);
                    printf(",%s", qdot_str[n][x]);
                }
                for(int x=0; x<3; x++) {
                    quadmath_snprintf(pdot_str[n][x], 512, "%.35Qf", pdot[n][x]);
                    printf(",%s", pdot_str[n][x]);
		}
            }

#else /* not QUAD */
            for(int x=0; x<3; x++) {
                printf(",%.35lf", q[n][x]);
            }
            for(int x=0; x<3; x++) {
                printf(",%.35lf", p[n][x]);
            }

            if (debug) {
                for(int x=0; x<3; x++) {
                    printf(",%lf", qdot[n][x]);
                }
                for(int x=0; x<3; x++) {
                    printf(",%lf", pdot[n][x]);
                }
            }
#endif /* QUAD */
        }
        
        
 if (debug) {
#ifdef QUAD

for(int x=0; x<4; x++) {
		quadmath_snprintf(H_str[x], 512, "%.35Qf", H[x]);
		printf(",%s", H_str[x]);
}
quadmath_snprintf(TotH_str, 512, "%.35Qf", H[0] + H[1] + H[2] + H[3]);
printf(",%s", TotH_str);

#else /* not QUAD */
	  for(int x=0; x<4; x++) {
                    printf(",%lf", H[x]);
}
printf(",%lf", H[0] + H[1] + H[2] + H[3]);
#endif /* QUAD */
}

        printf("\n");

        // Check if the maximum number of iterations will be exceeded.
        if((iteration >= iterations) && iterations !=0) {
            fprintf(stderr, "%s: exceeded maximum number of iterations \
                    (set to %i)\n", name, iterations);
            break;
        }
    }

    if (verbose) {
        // Record the end wall-time and calculate the difference.
#ifndef DARWIN
        clock_gettime(CLOCK_MONOTONIC_RAW, &end);
        uint64_t wall_time = ((end.tv_sec - start.tv_sec) * 1000 \
                + (end.tv_nsec - start.tv_nsec) / 1000000);

        fprintf(stderr,"%s: %" PRIu64 " milliseconds (wall-time)\n", name, wall_time);
#endif
    }

    // Re-run to check for more input.
    main(argc, argv);
    /************************************************************
     * END main
     ************************************************************/
}

void HamiltonEquations(PRECISION *H, PRECISION courant, bool calculate_hamiltonian) {
    /************************************************************
     * HamiltonEquations
     ************************************************************/
    /* This function calculates the right hand side of Hamilton's
     * Equations for Post-Minkowskian Gravity. This includes the
     * negative sign that appears in the equation for pdot.
     ************************************************************/

    // H = [Newtonian, H1, H2, H3]
    if (calculate_hamiltonian) H[0] = H[1] = H[2] = H[3] = ZERO;
    
    /************************************************************
     * COMPUTE: ps, En
     ************************************************************/
    for(int n=0; n<Ntot; n++) {
        ps[n] = P[n][0]*P[n][0] + P[n][1]*P[n][1] + P[n][2]*P[n][2];
        En[n] = sqrt(mass[n]*mass[n] + ps[n]);

        if(En[n] == 0) {
            fprintf(stderr, "%s: zero energy for particle %i.\n", name, n+1);
            exit(EXIT_FAILURE);
        }
    }

    /************************************************************
     * COMPUTE: r, Th, y, Xi
     ************************************************************/
    // Start with the assumption that the first two particles are closest.
    PRECISION rsmallest;
    // ats & bts are the labels for the closest particles.
    int ats=0, bts=1;

    for (int a=0; a<N; a++) {
        for (int b=0; b<N; b++) {
            if (b != a) {
                r[a][b] = sqrt((Q[b][0]-Q[a][0])*(Q[b][0]-Q[a][0]) 
                        + (Q[b][1]-Q[a][1])*(Q[b][1]-Q[a][1]) 
                        + (Q[b][2]-Q[a][2])*(Q[b][2]-Q[a][2]));

                if(r[a][b] == 0) {
                    fprintf(stderr, "%s: zero separation between particles %i & %i.\n", name, a+1, b+1);
                    exit(EXIT_FAILURE);
                }

                // This sets the value of rsmallest in the first iteration of the loop.
                rsmallest=r[ats][bts];

                // Check if the current pair of particles are closer.
                if(courant!=ZERO && (r[a][b] < rsmallest)) {
                    rsmallest = r[a][b];
                    ats = a;
                    bts = b;
                }

                Th[b][a] = (NEGATIVEONE/r[a][b]) * (P[b][0] * (Q[b][0]-Q[a][0]) \
                        + P[b][1] * (Q[b][1]-Q[a][1]) \
                        + P[b][2] * (Q[b][2]-Q[a][2]));
                y[b][a]  = (ONE/En[b]) * sqrt(mass[b]*mass[b] + Th[b][a]*Th[b][a]);
                Xi[a][b] = P[a][0]*P[b][0] \
                           + P[a][1]*P[b][1] \
                           + P[a][2]*P[b][2];
            }
        }
    }
    
    // Loops for test particles
    for (int a=0; a<N; a++) {
        for (int b=N; b<Ntot; b++) {
            if (b != a) {
                r[a][b] = sqrt((Q[b][0]-Q[a][0])*(Q[b][0]-Q[a][0]) 
                        + (Q[b][1]-Q[a][1])*(Q[b][1]-Q[a][1]) 
                        + (Q[b][2]-Q[a][2])*(Q[b][2]-Q[a][2]));
                
                r[b][a] = sqrt((Q[a][0]-Q[b][0])*(Q[a][0]-Q[b][0]) 
                        + (Q[a][1]-Q[b][1])*(Q[a][1]-Q[b][1]) 
                        + (Q[a][2]-Q[b][2])*(Q[a][2]-Q[b][2]));

                if(r[a][b] == 0) {
                    fprintf(stderr, "%s: zero separation between particles %i & %i.\n", name, a+1, b+1);
                    exit(EXIT_FAILURE);
                }
                
                if(r[b][a] == 0) {
                    fprintf(stderr, "%s: zero separation between particles %i & %i.\n", name, a+1, b+1);
                    exit(EXIT_FAILURE);
                }

                // This looks for the closest interacting particle/test particle pair
                // This implements adaptive time stepping for test particles
                if(courant<ZERO) {
                    rsmallest=r[ats][bts];

                    if(r[a][b] < rsmallest) {
                        rsmallest = r[a][b];
                        ats = a;
                        bts = b;
                    }
                
                    if(r[b][a] < rsmallest) {
                        rsmallest = r[b][a];
                        ats = a;
                        bts = b;
                    }
			    }

                Th[b][a] = (NEGATIVEONE/r[a][b]) * (P[b][0] * (Q[b][0]-Q[a][0]) \
                        + P[b][1] * (Q[b][1]-Q[a][1]) \
                        + P[b][2] * (Q[b][2]-Q[a][2]));
                y[b][a]  = (ONE/En[b]) * sqrt(mass[b]*mass[b] + Th[b][a]*Th[b][a]);
                Xi[a][b] = P[a][0]*P[b][0] \
                           + P[a][1]*P[b][1] \
                           + P[a][2]*P[b][2];
                           
                Th[a][b] = (NEGATIVEONE/r[b][a]) * (P[a][0] * (Q[a][0]-Q[b][0]) \
                        + P[a][1] * (Q[a][1]-Q[b][1]) \
                        + P[a][2] * (Q[a][2]-Q[b][2]));
                y[a][b]  = (ONE/En[a]) * sqrt(mass[a]*mass[a] + Th[a][b]*Th[a][b]);
                Xi[b][a] = P[b][0]*P[a][0] \
                           + P[b][1]*P[a][1] \
                           + P[b][2]*P[a][2];
            }
        }
    }

    /************************************************************
     * COMPUTE HAMILTON EQUATIONS
     ************************************************************/
    for(int n=0; n<Ntot; n++) {
        for (int x=0; x<3; x++) {
            pdot[n][x] = NEGATIVEONE * DHamiltonian(x,n);
            qdot[n][x] =             DHamiltonian(x+3,n);
        }
    }

    /************************************************************
     * CALCULATE TIMESTEP
     ************************************************************/
    if (courant!=ZERO) {
        PRECISION qdotsmallest = sqrt((qdot[bts][0]-qdot[ats][0]) * (qdot[bts][0]-qdot[ats][0]) \
                + (qdot[bts][1]-qdot[ats][1]) * (qdot[bts][1]-qdot[ats][1]) \
                + (qdot[bts][2]-qdot[ats][2]) * (qdot[bts][2]-qdot[ats][2]));

        PRECISION adapted_step = courant * (rsmallest/qdotsmallest);

        // Ensure the timestep is nonzero and less than the maximum size.
        if (adapted_step > ZERO){ 
            if(adapted_step < maxtimestep){
                timestep = adapted_step;
            }
            else{
                timestep=maxtimestep;
            }
        }
    }

    /************************************************************
     * HAMILTONIAN CALCULATOR
     ************************************************************/
    // Does not include contributions from test particles.
     
    if(calculate_hamiltonian) {

        // Subtract the rest mass from H1 to magnify changes in H.
        for(int n=0; n<N; n++) {
            H[1] += En[n] - mass[n];

            if(mass[n] > ZERO) {
                H[0] += HALF*ps[n]/mass[n];
            }
        }

        for(int a=0; a<N; a++) {
            for(int b=0; b<N; b++) {
                if(b != a) {
                    H[1] -= (gravity*HALF)*(En[a]*En[b]/r[a][b])*(ONE + (ps[a]/(En[a]*En[a])) + (ps[b]/(En[b]*En[b])));
                    H[2] += (gravity*QUARTER)*(((NEGATIVEONE)*(Th[a][b]*Th[b][a]) + SEVEN*Xi[a][b])/r[a][b]);

                    if(mass[b] > ZERO) {
                        H[3] -= (gravity*QUARTER)*(ONE/r[a][b])*(ONE/(En[a]*En[b]*(y[b][a]+ONE)*(y[b][a]+ONE)*y[b][a]))*( TWO*( TWO*Xi[a][b]*Xi[a][b]*Th[b][a]*Th[b][a] - TWO*(Th[a][b])*(-Th[b][a])*(Xi[a][b])*ps[b] + Th[a][b]*Th[a][b]*ps[b]*ps[b] - Xi[a][b]*Xi[a][b]*ps[b] )*(ONE/(En[b]*En[b])) + TWO*(-ps[a]*Th[b][a]*Th[b][a] + Th[a][b]*Th[a][b]*Th[b][a]*Th[b][a] + TWO*(Th[a][b])*(-Th[b][a])*(Xi[a][b]) + Xi[a][b]*Xi[a][b] - Th[a][b]*Th[a][b]*ps[b]) + ((-THREE*ps[a]*Th[b][a]*Th[b][a] + Th[a][b]*Th[a][b]*Th[b][a]*Th[b][a] + EIGHT*(Th[a][b])*(-Th[b][a])*(Xi[a][b]) + ps[a]*ps[b] - THREE*Th[a][b]*Th[a][b]*ps[b] )*y[b][a] ) );
                        H[0] -= HALF*gravity*mass[a]*mass[b]/r[a][b];
                    }
                    else {
                        H[3] -= (gravity*QUARTER)*(ps[b]*Th[a][b]*Th[a][b]*(NEGATIVEONE + y[b][a])*(THREE + y[b][a]) - ps[a]*ps[b]*(ONE + y[b][a])*(NEGATIVEONE + THREE*y[b][a]) + FOUR*Xi[a][b]*(-TWO*Th[a][b]*Th[b][a] + Xi[a][b]*y[b][a]))/(sqrt(ps[a])*sqrt(ps[b])*r[a][b]*(ONE + y[b][a])*(ONE + y[b][a]));
                    }
                }
            }
        }

        switch (hamiltonian_to_use) {
            case 0:        H[1]=H[2]=H[3]=ZERO; break;
            case 1:   H[0]     =H[2]=H[3]=ZERO; break;
            case 12:  H[0]          =H[3]=ZERO; break;
            case 13:  H[0]     =H[2]     =ZERO; break;
            default: case 123: H[0]               =ZERO; break;
        }
    }

}

PRECISION DHamiltonian(unsigned int zindex, unsigned const int c)
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
     * dps, dEn, dr, dy, dTh, dXi are derivatives of ps, En, m, r, y, Th, Xi. 
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
		dps[c] = ZERO;
        dEn[c] = ZERO;
		}
		
        for(int a=0; a<N; a++) {
            dps[a] = ZERO;
            dEn[a] = ZERO;

            if(a!=c) {
                dr[a][c]  = (ONE/r[a][c]) * (Q[c][zindex]-Q[a][zindex]);
                dTh[a][c] = (ONE/r[a][c]) * (P[a][zindex]-Th[a][c]*dr[a][c]);
                dXi[a][c] = ZERO;

                if(mass[a] > ZERO) {
                    dy[a][c] = (En[a]*dTh[a][c]*Th[a][c] + dEn[a]*(-mass[a]*mass[a] - Th[a][c]*Th[a][c]))/(En[a]*En[a]*sqrt(mass[a]*mass[a] + Th[a][c]*Th[a][c]));
                }
                else {
                    dy[a][c] = Sign(Th[a][c])*(En[a]*dTh[a][c] - dEn[a]*Th[a][c])/(En[a]*En[a]);
                }
            }
        }
        for(int b=0; b<N; b++) {
            if(b!=c) {
                dr[c][b]  = (ONE/r[c][b]) * (Q[b][zindex]-Q[c][zindex]) * (NEGATIVEONE);
                dTh[c][b] = (ONE/r[c][b]) * (P[c][zindex]*(NEGATIVEONE) - Th[c][b]*dr[c][b]);
                dXi[c][b] = ZERO;

                if(mass[c] > ZERO) {
                    dy[c][b] = (En[c]*dTh[c][b]*Th[c][b] + dEn[c]*(-mass[c]*mass[c] - Th[c][b]*Th[c][b]))/(En[c]*En[c]*sqrt(mass[c]*mass[c] + Th[c][b]*Th[c][b]));
                }
                else {
                    dy[c][b] = Sign(Th[c][b])*(En[c]*dTh[c][b] - dEn[c]*Th[c][b])/(En[c]*En[c]);
                }
            }
        }
    }

    else if(zindex > 2) {
        // This allows us to use zindex directly as a fixed subscript.
        zindex -= 3;

        for(int a=0; a<N; a++) {
            if(a==c) {
                dps[a] = TWO*P[a][zindex];
            }
            else {
                dps[a] = ZERO;
            }

            dEn[a] = (HALF/(En[a])) * dps[a];
        }
        
            // For Test particles
            if(c>=N) {
			dps[c] = TWO*P[c][zindex];
            dEn[c] = (HALF/(En[c])) * dps[c];
			}

        for(int a=0; a<N; a++) {
            if(a!=c) {
                dr[a][c]  = ZERO;
                dTh[a][c] = ZERO;
                dXi[a][c] = P[a][zindex];

                if(mass[a] > ZERO) {
                    dy[a][c] = (En[a]*dTh[a][c]*Th[a][c] + dEn[a]*(-mass[a]*mass[a] - Th[a][c]*Th[a][c]))/(En[a]*En[a]*sqrt(mass[a]*mass[a] + Th[a][c]*Th[a][c]));
                }
                else {
                    dy[a][c] = Sign(Th[a][c])*(En[a]*dTh[a][c] - dEn[a]*Th[a][c])/(En[a]*En[a]);
                }
            }
        }

        for(int b=0; b<N; b++) {
            if(b!=c) {
                dr[c][b]  = ZERO;
                dTh[c][b] = (ONE/r[c][b]) * (Q[b][zindex]-Q[c][zindex]);
                dXi[c][b] = P[b][zindex];

                if(mass[c] > ZERO) {
                    dy[c][b] = (En[c]*dTh[c][b]*Th[c][b] + dEn[c]*(-mass[c]*mass[c] - Th[c][b]*Th[c][b]))/(En[c]*En[c]*sqrt(mass[c]*mass[c] + Th[c][b]*Th[c][b]));
                }
                else {
                    dy[c][b] = Sign(Th[c][b])*(En[c]*dTh[c][b] - dEn[c]*Th[c][b])/(En[c]*En[c]);
                }
            }
        }
    }

    /************************************************************
     * CALCULATE dH (the derivative of the hamiltonian)
     ************************************************************/
    // dH[4] = {dHnewt,dH1,dH2,dH3}
    PRECISION dH[4] = {ZERO,ZERO,ZERO,ZERO};

    for(int n=0; n<N; n++)  {
        dH[1] += dEn[n];
        dH[0] += HALF*dps[n]/mass[n];
        
        // For Test particles
        if(c>=N) {
		dH[1] += dEn[c];
        dH[0] += HALF*dps[c]/mass[c];
		}
    }

    for(int a=0; a<N; a++) {
        if(c!=a) {
            dH[1] -= (gravity*HALF)*((-(En[a]*En[c]*(En[c]*En[c]*ps[a] + En[a]*En[a]*(En[c]*En[c] + ps[c]))*dr[a][c]) + (-(dEn[a]*En[c]*En[c]*En[c]*ps[a]) + En[a]*En[c]*En[c]*(dps[a]*En[c] + dEn[c]*ps[a]) + En[a]*En[a]*En[a]*(dps[c]*En[c] + dEn[c]*(En[c]*En[c] - ps[c])) + dEn[a]*En[a]*En[a]*En[c]*(En[c]*En[c] + ps[c]))*r[a][c])/(En[a]*En[a]*En[c]*En[c]*r[a][c]*r[a][c]));
            dH[2] += (gravity*QUARTER)*((SEVEN*dXi[a][c]*r[a][c] - dTh[c][a]*r[a][c]*Th[a][c] - dTh[a][c]*r[a][c]*Th[c][a] + dr[a][c]*Th[a][c]*Th[c][a] - SEVEN*dr[a][c]*Xi[a][c])/(r[a][c]*r[a][c]));

            if(mass[c] > ZERO) {
                dH[3] -= (gravity*QUARTER)*(-((TWO*En[a]*En[c]*dy[c][a]*r[a][c]*y[c][a]*(TWO*ps[c]*ps[c]*Th[a][c]*Th[a][c] + FOUR*ps[c]*Th[a][c]*Th[c][a]*Xi[a][c] - TWO*(ps[c] - TWO*Th[c][a]*Th[c][a])*Xi[a][c]*Xi[a][c] + En[c]*En[c]*(TWO*(-(Th[a][c]*Th[c][a]) + Xi[a][c])*(-(Th[a][c]*Th[c][a]) + Xi[a][c]) + Th[a][c]*Th[c][a]*(Th[a][c]*Th[c][a] - EIGHT*Xi[a][c])*y[c][a] - ps[a]*Th[c][a]*Th[c][a]*(TWO + THREE*y[c][a]) + ps[c]*(ps[a]*y[c][a] - Th[a][c]*Th[a][c]*(TWO + THREE*y[c][a])))) + dEn[c]*En[a]*r[a][c]*y[c][a]*(ONE + y[c][a])*(TWO*ps[c]*ps[c]*Th[a][c]*Th[a][c] + FOUR*ps[c]*Th[a][c]*Th[c][a]*Xi[a][c] - TWO*(ps[c] - TWO*Th[c][a]*Th[c][a])*Xi[a][c]*Xi[a][c] + En[c]*En[c]*(TWO*(-(Th[a][c]*Th[c][a]) + Xi[a][c])*(-(Th[a][c]*Th[c][a]) + Xi[a][c]) + Th[a][c]*Th[c][a]*(Th[a][c]*Th[c][a] - EIGHT*Xi[a][c])*y[c][a] - ps[a]*Th[c][a]*Th[c][a]*(TWO + THREE*y[c][a]) + ps[c]*(ps[a]*y[c][a] - Th[a][c]*Th[a][c]*(TWO + THREE*y[c][a])))) - En[a]*En[c]*dy[c][a]*r[a][c]*(ONE + y[c][a])*(-TWO*ps[c]*ps[c]*Th[a][c]*Th[a][c] - FOUR*Th[c][a]*Th[c][a]*Xi[a][c]*Xi[a][c] + TWO*ps[c]*Xi[a][c]*(-TWO*Th[a][c]*Th[c][a] + Xi[a][c]) + En[c]*En[c]*(-TWO*(-(Th[a][c]*Th[c][a]) + Xi[a][c])*(-(Th[a][c]*Th[c][a]) + Xi[a][c]) + Th[a][c]*Th[c][a]*(-(Th[a][c]*Th[c][a]) + EIGHT*Xi[a][c])*y[c][a] + ps[a]*Th[c][a]*Th[c][a]*(TWO + THREE*y[c][a]) + ps[c]*(-(ps[a]*y[c][a]) + Th[a][c]*Th[a][c]*(TWO + THREE*y[c][a])))) - En[a]*En[c]*dr[a][c]*y[c][a]*(ONE + y[c][a])*(-TWO*ps[c]*ps[c]*Th[a][c]*Th[a][c] - FOUR*Th[c][a]*Th[c][a]*Xi[a][c]*Xi[a][c] + TWO*ps[c]*Xi[a][c]*(-TWO*Th[a][c]*Th[c][a] + Xi[a][c]) + En[c]*En[c]*(-TWO*(-(Th[a][c]*Th[c][a]) + Xi[a][c])*(-(Th[a][c]*Th[c][a]) + Xi[a][c]) + Th[a][c]*Th[c][a]*(-(Th[a][c]*Th[c][a]) + EIGHT*Xi[a][c])*y[c][a] + ps[a]*Th[c][a]*Th[c][a]*(TWO + THREE*y[c][a]) + ps[c]*(-(ps[a]*y[c][a]) + Th[a][c]*Th[a][c]*(TWO + THREE*y[c][a])))) - dEn[a]*En[c]*r[a][c]*y[c][a]*(ONE + y[c][a])*(-TWO*ps[c]*ps[c]*Th[a][c]*Th[a][c] - FOUR*Th[c][a]*Th[c][a]*Xi[a][c]*Xi[a][c] + TWO*ps[c]*Xi[a][c]*(-TWO*Th[a][c]*Th[c][a] + Xi[a][c]) + En[c]*En[c]*(-TWO*(-(Th[a][c]*Th[c][a]) + Xi[a][c])*(-(Th[a][c]*Th[c][a]) + Xi[a][c]) + Th[a][c]*Th[c][a]*(-(Th[a][c]*Th[c][a]) + EIGHT*Xi[a][c])*y[c][a] + ps[a]*Th[c][a]*Th[c][a]*(TWO + THREE*y[c][a]) + ps[c]*(-(ps[a]*y[c][a]) + Th[a][c]*Th[a][c]*(TWO + THREE*y[c][a])))) + En[a]*r[a][c]*y[c][a]*(ONE + y[c][a])*(TWO*dps[c]*En[c]*En[c]*En[c]*Th[a][c]*Th[a][c] + FOUR*ps[c]*ps[c]*Th[a][c]*(-(En[c]*dTh[a][c]) + dEn[c]*Th[a][c]) - FOUR*En[c]*En[c]*En[c]*dXi[a][c]*Xi[a][c] + FOUR*En[c]*En[c]*En[c]*dTh[c][a]*Th[a][c]*Xi[a][c] + TWO*dps[c]*En[c]*Xi[a][c]*Xi[a][c] - dps[c]*En[c]*En[c]*En[c]*ps[a]*y[c][a] + THREE*dps[c]*En[c]*En[c]*En[c]*Th[a][c]*Th[a][c]*y[c][a] + EIGHT*En[c]*En[c]*En[c]*dTh[c][a]*Th[a][c]*Xi[a][c]*y[c][a] + Th[c][a]*Th[c][a]*(-EIGHT*En[c]*dXi[a][c]*Xi[a][c] + EIGHT*dEn[c]*Xi[a][c]*Xi[a][c] + En[c]*En[c]*En[c]*(dy[c][a]*(THREE*ps[a] - Th[a][c]*Th[a][c]) - TWO*dTh[a][c]*Th[a][c]*(TWO + y[c][a]) + dps[a]*(TWO + THREE*y[c][a]))) + ps[c]*((-FOUR*dps[c]*En[c] + THREE*En[c]*En[c]*En[c]*dy[c][a])*Th[a][c]*Th[a][c] + FOUR*En[c]*(dXi[a][c] - dTh[a][c]*Th[c][a])*Xi[a][c] - FOUR*dEn[c]*Xi[a][c]*Xi[a][c] - En[c]*En[c]*En[c]*(ps[a]*dy[c][a] + dps[a]*y[c][a]) + Th[a][c]*(-FOUR*En[c]*dTh[c][a]*Xi[a][c] + Th[c][a]*(-FOUR*En[c]*dXi[a][c] + EIGHT*dEn[c]*Xi[a][c]) + TWO*En[c]*En[c]*En[c]*dTh[a][c]*(TWO + THREE*y[c][a]))) + Th[c][a]*(-FOUR*En[c]*Xi[a][c]*(dps[c]*Th[a][c] + TWO*dTh[c][a]*Xi[a][c]) + En[c]*En[c]*En[c]*(FOUR*dXi[a][c]*Th[a][c]*(ONE + TWO*y[c][a]) + FOUR*Xi[a][c]*(TWO*dy[c][a]*Th[a][c] + dTh[a][c]*(ONE + TWO*y[c][a])) + TWO*dTh[c][a]*(-(Th[a][c]*Th[a][c]*(TWO + y[c][a])) + ps[a]*(TWO + THREE*y[c][a]))))))/(En[a]*En[a]*En[c]*En[c]*En[c]*En[c]*r[a][c]*r[a][c]*y[c][a]*y[c][a]*(ONE + y[c][a])*(ONE + y[c][a])*(ONE + y[c][a]))));
                dH[0] += HALF*gravity*mass[a]*mass[c]*dr[a][c]/(r[a][c]*r[a][c]);
            }
            else {
                PRECISION signTh = Sign(Th[c][a]);
                dH[3] -= (gravity*QUARTER)*(TWO*En[a]*ps[c]*dr[a][c]*(ONE + y[c][a])*(EIGHT*signTh*sqrt(ps[c])*Th[a][c]*Xi[a][c]*y[c][a] - FOUR*Xi[a][c]*Xi[a][c]*y[c][a] - ps[c]*Th[a][c]*Th[a][c]*(NEGATIVEONE + y[c][a])*(THREE + y[c][a]) +  ps[a]*ps[c]*(ONE + y[c][a])*(NEGATIVEONE + THREE*y[c][a])) + r[a][c]*(dps[c]*En[a]*(ONE + y[c][a])*(EIGHT*signTh*sqrt(ps[c])*Th[a][c]*Xi[a][c]*y[c][a] - TWO*SIX*Xi[a][c]*Xi[a][c]*y[c][a] - ps[c]*Th[a][c]*Th[a][c]*(THREE + y[c][a]*(TWO + y[c][a])) + ps[a]*ps[c]*(ONE + y[c][a]*(TWO + THREE*y[c][a]))) - TWO*(dps[a]*En[a]*ps[c]*ps[c]*(ONE + y[c][a])*(ONE + y[c][a])*(NEGATIVEONE + THREE*y[c][a]) + dEn[a]*ps[c]*(ONE + y[c][a])*(-EIGHT*signTh*sqrt(ps[c])*Th[a][c]*Xi[a][c]*y[c][a] + FOUR*Xi[a][c]*Xi[a][c]*y[c][a] + ps[c]*Th[a][c]*Th[a][c]*(NEGATIVEONE + y[c][a])*(THREE + y[c][a]) - ps[a]*ps[c]*(ONE + y[c][a])*(NEGATIVEONE + THREE*y[c][a])) + TWO*En[a]*sqrt(ps[c])*(dTh[c][a]*(ONE + y[c][a])*(FOUR*sqrt(ps[c])*Th[a][c]*Xi[a][c] - FOUR*signTh*Xi[a][c]*Xi[a][c] - signTh*ps[c]*Th[a][c]*Th[a][c]*(TWO + y[c][a]) + signTh*ps[a]*ps[c]*(TWO + THREE*y[c][a])) + sqrt(ps[c])*(FOUR*dXi[a][c]*(signTh*sqrt(ps[c])*Th[a][c] - Xi[a][c])*y[c][a]*(ONE + y[c][a]) - sqrt(ps[c])*dTh[a][c]*(ONE + y[c][a])*(-FOUR*signTh*Xi[a][c]*y[c][a] + sqrt(ps[c])*Th[a][c]*(NEGATIVEONE + y[c][a])*(THREE + y[c][a])) + dy[c][a]*(-EIGHT*signTh*sqrt(ps[c])*Th[a][c]*Xi[a][c]*y[c][a] + TWO*Xi[a][c]*Xi[a][c]*(ONE + THREE*y[c][a]) + ps[c]*(-THREE*ps[a]*y[c][a]*(ONE + y[c][a]) + Th[a][c]*Th[a][c]*(-TWO + y[c][a]*(THREE + y[c][a])))))))))/(TWO*En[a]*En[a]*sqrt(ps[c])*ps[c]*r[a][c]*r[a][c]*(ONE + y[c][a])*(ONE + y[c][a])*(ONE + y[c][a]));
            }
        }
    }
    for(int b=0; b<N; b++) {
        if(b!=c) {
            dH[1] -= (gravity*HALF)*((-(En[c]*En[b]*(En[b]*En[b]*ps[c] + En[c]*En[c]*(En[b]*En[b] + ps[b]))*dr[c][b]) + (-(dEn[c]*En[b]*En[b]*En[b]*ps[c]) + En[c]*En[b]*En[b]*(dps[c]*En[b] + dEn[b]*ps[c]) + En[c]*En[c]*En[c]*(dps[b]*En[b] + dEn[b]*(En[b]*En[b] - ps[b])) + dEn[c]*En[c]*En[c]*En[b]*(En[b]*En[b] + ps[b]))*r[c][b])/(En[c]*En[c]*En[b]*En[b]*r[c][b]*r[c][b]));
            dH[2] += (gravity*QUARTER)*((SEVEN*dXi[c][b]*r[c][b] - dTh[b][c]*r[c][b]*Th[c][b] - dTh[c][b]*r[c][b]*Th[b][c] + dr[c][b]*Th[c][b]*Th[b][c] - SEVEN*dr[c][b]*Xi[c][b])/(r[c][b]*r[c][b]));

            if(mass[b] > ZERO) {
                dH[3] -= (gravity*QUARTER)*(-((TWO*En[c]*En[b]*dy[b][c]*r[c][b]*y[b][c]*(TWO*ps[b]*ps[b]*Th[c][b]*Th[c][b] + FOUR*ps[b]*Th[c][b]*Th[b][c]*Xi[c][b] - TWO*(ps[b] - TWO*Th[b][c]*Th[b][c])*Xi[c][b]*Xi[c][b] + En[b]*En[b]*(TWO*(-(Th[c][b]*Th[b][c]) + Xi[c][b])*(-(Th[c][b]*Th[b][c]) + Xi[c][b]) + Th[c][b]*Th[b][c]*(Th[c][b]*Th[b][c] - EIGHT*Xi[c][b])*y[b][c] - ps[c]*Th[b][c]*Th[b][c]*(TWO + THREE*y[b][c]) + ps[b]*(ps[c]*y[b][c] - Th[c][b]*Th[c][b]*(TWO + THREE*y[b][c])))) + dEn[b]*En[c]*r[c][b]*y[b][c]*(ONE + y[b][c])*(TWO*ps[b]*ps[b]*Th[c][b]*Th[c][b] + FOUR*ps[b]*Th[c][b]*Th[b][c]*Xi[c][b] - TWO*(ps[b] - TWO*Th[b][c]*Th[b][c])*Xi[c][b]*Xi[c][b] + En[b]*En[b]*(TWO*(-(Th[c][b]*Th[b][c]) + Xi[c][b])*(-(Th[c][b]*Th[b][c]) + Xi[c][b]) + Th[c][b]*Th[b][c]*(Th[c][b]*Th[b][c] - EIGHT*Xi[c][b])*y[b][c] - ps[c]*Th[b][c]*Th[b][c]*(TWO + THREE*y[b][c]) + ps[b]*(ps[c]*y[b][c] - Th[c][b]*Th[c][b]*(TWO + THREE*y[b][c])))) - En[c]*En[b]*dy[b][c]*r[c][b]*(ONE + y[b][c])*(-TWO*ps[b]*ps[b]*Th[c][b]*Th[c][b] - FOUR*Th[b][c]*Th[b][c]*Xi[c][b]*Xi[c][b] + TWO*ps[b]*Xi[c][b]*(-TWO*Th[c][b]*Th[b][c] + Xi[c][b]) + En[b]*En[b]*(-TWO*(-(Th[c][b]*Th[b][c]) + Xi[c][b])*(-(Th[c][b]*Th[b][c]) + Xi[c][b]) + Th[c][b]*Th[b][c]*(-(Th[c][b]*Th[b][c]) + EIGHT*Xi[c][b])*y[b][c] + ps[c]*Th[b][c]*Th[b][c]*(TWO + THREE*y[b][c]) + ps[b]*(-(ps[c]*y[b][c]) + Th[c][b]*Th[c][b]*(TWO + THREE*y[b][c])))) - En[c]*En[b]*dr[c][b]*y[b][c]*(ONE + y[b][c])*(-TWO*ps[b]*ps[b]*Th[c][b]*Th[c][b] - FOUR*Th[b][c]*Th[b][c]*Xi[c][b]*Xi[c][b] + TWO*ps[b]*Xi[c][b]*(-TWO*Th[c][b]*Th[b][c] + Xi[c][b]) + En[b]*En[b]*(-TWO*(-(Th[c][b]*Th[b][c]) + Xi[c][b])*(-(Th[c][b]*Th[b][c]) + Xi[c][b]) + Th[c][b]*Th[b][c]*(-(Th[c][b]*Th[b][c]) + EIGHT*Xi[c][b])*y[b][c] + ps[c]*Th[b][c]*Th[b][c]*(TWO + THREE*y[b][c]) + ps[b]*(-(ps[c]*y[b][c]) + Th[c][b]*Th[c][b]*(TWO + THREE*y[b][c])))) - dEn[c]*En[b]*r[c][b]*y[b][c]*(ONE + y[b][c])*(-TWO*ps[b]*ps[b]*Th[c][b]*Th[c][b] - FOUR*Th[b][c]*Th[b][c]*Xi[c][b]*Xi[c][b] + TWO*ps[b]*Xi[c][b]*(-TWO*Th[c][b]*Th[b][c] + Xi[c][b]) + En[b]*En[b]*(-TWO*(-(Th[c][b]*Th[b][c]) + Xi[c][b])*(-(Th[c][b]*Th[b][c]) + Xi[c][b]) + Th[c][b]*Th[b][c]*(-(Th[c][b]*Th[b][c]) + EIGHT*Xi[c][b])*y[b][c] + ps[c]*Th[b][c]*Th[b][c]*(TWO + THREE*y[b][c]) + ps[b]*(-(ps[c]*y[b][c]) + Th[c][b]*Th[c][b]*(TWO + THREE*y[b][c])))) + En[c]*r[c][b]*y[b][c]*(ONE + y[b][c])*(TWO*dps[b]*En[b]*En[b]*En[b]*Th[c][b]*Th[c][b] + FOUR*ps[b]*ps[b]*Th[c][b]*(-(En[b]*dTh[c][b]) + dEn[b]*Th[c][b]) - FOUR*En[b]*En[b]*En[b]*dXi[c][b]*Xi[c][b] + FOUR*En[b]*En[b]*En[b]*dTh[b][c]*Th[c][b]*Xi[c][b] + TWO*dps[b]*En[b]*Xi[c][b]*Xi[c][b] - dps[b]*En[b]*En[b]*En[b]*ps[c]*y[b][c] + THREE*dps[b]*En[b]*En[b]*En[b]*Th[c][b]*Th[c][b]*y[b][c] + EIGHT*En[b]*En[b]*En[b]*dTh[b][c]*Th[c][b]*Xi[c][b]*y[b][c] + Th[b][c]*Th[b][c]*(-EIGHT*En[b]*dXi[c][b]*Xi[c][b] + EIGHT*dEn[b]*Xi[c][b]*Xi[c][b] + En[b]*En[b]*En[b]*(dy[b][c]*(THREE*ps[c] - Th[c][b]*Th[c][b]) - TWO*dTh[c][b]*Th[c][b]*(TWO + y[b][c]) + dps[c]*(TWO + THREE*y[b][c]))) + ps[b]*((-FOUR*dps[b]*En[b] + THREE*En[b]*En[b]*En[b]*dy[b][c])*Th[c][b]*Th[c][b] + FOUR*En[b]*(dXi[c][b] - dTh[c][b]*Th[b][c])*Xi[c][b] - FOUR*dEn[b]*Xi[c][b]*Xi[c][b] - En[b]*En[b]*En[b]*(ps[c]*dy[b][c] + dps[c]*y[b][c]) + Th[c][b]*(-FOUR*En[b]*dTh[b][c]*Xi[c][b] + Th[b][c]*(-FOUR*En[b]*dXi[c][b] + EIGHT*dEn[b]*Xi[c][b]) + TWO*En[b]*En[b]*En[b]*dTh[c][b]*(TWO + THREE*y[b][c]))) + Th[b][c]*(-FOUR*En[b]*Xi[c][b]*(dps[b]*Th[c][b] + TWO*dTh[b][c]*Xi[c][b]) + En[b]*En[b]*En[b]*(FOUR*dXi[c][b]*Th[c][b]*(ONE + TWO*y[b][c]) + FOUR*Xi[c][b]*(TWO*dy[b][c]*Th[c][b] + dTh[c][b]*(ONE + TWO*y[b][c])) + TWO*dTh[b][c]*(-(Th[c][b]*Th[c][b]*(TWO + y[b][c])) + ps[c]*(TWO + THREE*y[b][c]))))))/(En[c]*En[c]*En[b]*En[b]*En[b]*En[b]*r[c][b]*r[c][b]*y[b][c]*y[b][c]*(ONE + y[b][c])*(ONE + y[b][c])*(ONE + y[b][c]))));
                dH[0] += HALF*gravity*mass[c]*mass[b]*dr[c][b]/(r[c][b]*r[c][b]);
            }
            else {
                PRECISION signTh = Sign(Th[b][c]);
                dH[3] -= (gravity*QUARTER)*(TWO*En[c]*ps[b]*dr[c][b]*(ONE + y[b][c])*(EIGHT*signTh*sqrt(ps[b])*Th[c][b]*Xi[c][b]*y[b][c] - FOUR*Xi[c][b]*Xi[c][b]*y[b][c] - ps[b]*Th[c][b]*Th[c][b]*(NEGATIVEONE + y[b][c])*(THREE + y[b][c]) +  ps[c]*ps[b]*(ONE + y[b][c])*(NEGATIVEONE + THREE*y[b][c])) + r[c][b]*(dps[b]*En[c]*(ONE + y[b][c])*(EIGHT*signTh*sqrt(ps[b])*Th[c][b]*Xi[c][b]*y[b][c] - TWO*SIX*Xi[c][b]*Xi[c][b]*y[b][c] - ps[b]*Th[c][b]*Th[c][b]*(THREE + y[b][c]*(TWO + y[b][c])) + ps[c]*ps[b]*(ONE + y[b][c]*(TWO + THREE*y[b][c]))) - TWO*(dps[c]*En[c]*ps[b]*ps[b]*(ONE + y[b][c])*(ONE + y[b][c])*(NEGATIVEONE + THREE*y[b][c]) + dEn[c]*ps[b]*(ONE + y[b][c])*(-EIGHT*signTh*sqrt(ps[b])*Th[c][b]*Xi[c][b]*y[b][c] + FOUR*Xi[c][b]*Xi[c][b]*y[b][c] + ps[b]*Th[c][b]*Th[c][b]*(NEGATIVEONE + y[b][c])*(THREE + y[b][c]) - ps[c]*ps[b]*(ONE + y[b][c])*(NEGATIVEONE + THREE*y[b][c])) + TWO*En[c]*sqrt(ps[b])*(dTh[b][c]*(ONE + y[b][c])*(FOUR*sqrt(ps[b])*Th[c][b]*Xi[c][b] - FOUR*signTh*Xi[c][b]*Xi[c][b] - signTh*ps[b]*Th[c][b]*Th[c][b]*(TWO + y[b][c]) + signTh*ps[c]*ps[b]*(TWO + THREE*y[b][c])) + sqrt(ps[b])*(FOUR*dXi[c][b]*(signTh*sqrt(ps[b])*Th[c][b] - Xi[c][b])*y[b][c]*(ONE + y[b][c]) - sqrt(ps[b])*dTh[c][b]*(ONE + y[b][c])*(-FOUR*signTh*Xi[c][b]*y[b][c] + sqrt(ps[b])*Th[c][b]*(NEGATIVEONE + y[b][c])*(THREE + y[b][c])) + dy[b][c]*(-EIGHT*signTh*sqrt(ps[b])*Th[c][b]*Xi[c][b]*y[b][c] + TWO*Xi[c][b]*Xi[c][b]*(ONE + THREE*y[b][c]) + ps[b]*(-THREE*ps[c]*y[b][c]*(ONE + y[b][c]) + Th[c][b]*Th[c][b]*(-TWO + y[b][c]*(THREE + y[b][c])))))))))/(TWO*En[c]*En[c]*sqrt(ps[b])*ps[b]*r[c][b]*r[c][b]*(ONE + y[b][c])*(ONE + y[b][c])*(ONE + y[b][c]));
            }
        }
    }


    switch (hamiltonian_to_use) {
        case 0:         dH[1]=dH[2]=dH[3]=ZERO; break;
        case 1:   dH[0]      =dH[2]=dH[3]=ZERO; break;
        case 12:  dH[0]            =dH[3]=ZERO; break;
        case 13:  dH[0]      =dH[2]      =ZERO; break;
        default: case 123: dH[0]                  =ZERO; break;
    }

    return(dH[0]+dH[1]+dH[2]+dH[3]);
}

PRECISION Sign(PRECISION A) {
    /* Return the sign of A. */
    if(A == ZERO)
        return(ZERO);
    else  {
        if(signbit(A)==true) {
            return(NEGATIVEONE);
        }
        else {
            return(ONE);
        }
    }
}

/************************************************************
 * END pomin.c
 ************************************************************/
