# include "bioplib/pdb.h"
# include "arrays.h"

#ifndef REGRESSION 
#define REGRESSION

void draw_regression_line(double **coordinates,
                          double *eigenVector,
                          int numberOfPoints,
                          char *chainLabel,
                          FILE *wfp);


void compute_best_fit_line(double **coordinates,
			   int numberOfPoints,
			   int numberOfDimensions,
			   double *centroid,
			   double *eigenVector);

#endif
