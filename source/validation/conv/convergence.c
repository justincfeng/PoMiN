/************************************************************
 * convergence.c - calculate convergence values from the output of autorun.sh
 ************************************************************/
/* This file calculates the convergence values from a csv file in the
 * form of <output-file>.csv from pomin(1). It takes the
 * <output-file> csv file generated from autorun.sh as its single
 * argument. It is written in C to take advantage of the quadmath
 * precision.
 ************************************************************/

/************************************************************
 * HEADERS
 ************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <quadmath.h>

#define QX2 8
#define QY2 9
#define QZ2 10
#define PX2 11
#define PY2 12
#define PZ2 13

/************************************************************
 * MAIN
 ************************************************************/
int main(int argc, char *argv[]) {
  /* Throw an error if there is not just one command line argument. */
  if (argc != 2 || (argv[argc-1] && strcmp(argv[argc-1],"-h") == 0)) {
      /* argv[0] is the program name. */
      fprintf(stderr,"USAGE: %s <csv-file>\n", argv[0]);
      exit(EXIT_FAILURE);
  }

  /************************************************************
   * SETUP
   ************************************************************/
  /* Define 2d-array to store csv values as floats: q[row][column]. */
  __float128 q[512][512];
  /* Define 1d-array to store momenta squared */
  __float128 psquared[256];
  /* Define 1d-array to store final convergence values. */
  __float128 Q[256];
  /* Define variables to hold csv size. */
  int numrow=0, numcolm=0;

  /* Open input csv file. */
  FILE* infile = fopen(argv[1], "r");
  /* Exit with error if file not found. */
  if (infile == NULL) {
      fprintf(stderr,"%s: ERROR: Could not open csv file\n",argv[0]);
      exit(EXIT_FAILURE);
  }

  /* Define variable to hold rows to process. */
  char row[10000];              /* TODO: make row arbitrarily long */

  /* Determine number of csv columns. */
  if(fgets (row, sizeof(row), infile)!=NULL) {

      for(int i=0; row[i] !=0; i++) {
          if(row[i] == ',') {
              numcolm++;
          }
      }
  }

  /************************************************************
   * PULL DATA FROM CSV
   ************************************************************/
  /* Iterate over csv rows. */
  for(int r=0; fgets(row, sizeof(row), infile); r++) {

      /* Convert each entry from a string into a float. */
      q[r][0] = strtoflt128(strtok(row,","), NULL);

      /* Iterate over csv columns. */
      for(int c=1; c<=numcolm; c++) {
          q[r][c] = strtoflt128(strtok(NULL,","), NULL);
      }
      /* Catch the number of rows in a variable for later use. */
      numrow++;
  }
  fclose(infile);

  /************************************************************
   * PRINT RESULTS
   ************************************************************/
  /* Define a recycled buffer to print quadmath strings with. */
  char buf[256];

  /* Iterate back through the rows (of q this time). */
  for(int b=0;b<=numrow-3;b++) {
      /* Calculate the convergence value. */
      //Q[b] = sqrtq( powq((q[b+0][8] - q[b+1][8]),2) + powq((q[b+0][9] - q[b+1][9]),2) )
      //     / sqrtq( powq((q[b+1][8] - q[b+2][8]),2) + powq((q[b+1][9] - q[b+2][9]),2) );

      psquared[b+0] = q[b+0][PX2] * q[b+0][PX2] + q[b+0][PY2] * q[b+0][PY2] + q[b+0][PZ2] * q[b+0][PZ2];
      psquared[b+1] = q[b+1][PX2] * q[b+1][PX2] + q[b+1][PY2] * q[b+1][PY2] + q[b+1][PZ2] * q[b+1][PZ2];
      psquared[b+2] = q[b+2][PX2] * q[b+2][PX2] + q[b+2][PY2] * q[b+2][PY2] + q[b+2][PZ2] * q[b+2][PZ2];

      Q[b] = fabsq( ( psquared[b+0] - psquared[b+1] ) / ( psquared[b+1] - psquared[b+2] ) );

      /* Print this value. */
      quadmath_snprintf(buf, sizeof(buf), "%.35Qf", Q[b]);
      printf("%s\n", buf);
  }

  exit(EXIT_SUCCESS);
}

/************************************************************
 * END convergence.c
 ************************************************************/
