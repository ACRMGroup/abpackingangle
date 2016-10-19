#include "arrays.h"

#ifndef REGRESSION 
#define REGRESSION

void DrawRegressionLine(REAL **coordinates,
                          REAL *eigenVector,
                          int numberOfPoints,
                          char *chainLabel,
                          FILE *wfp);


void ComputeBestFitLine(REAL **coordinates,
			   int numberOfPoints,
			   int numberOfDimensions,
			   REAL *centroid,
			   REAL *eigenVector);

#endif
