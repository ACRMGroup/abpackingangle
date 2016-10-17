#include "arrays.h"

#ifndef REGRESSION 
#define REGRESSION

void draw_regression_line(REAL **coordinates,
                          REAL *eigenVector,
                          int numberOfPoints,
                          char *chainLabel,
                          FILE *wfp);


void compute_best_fit_line(REAL **coordinates,
			   int numberOfPoints,
			   int numberOfDimensions,
			   REAL *centroid,
			   REAL *eigenVector);

#endif
