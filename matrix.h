# include <stdio.h>
# include <string.h>
# include "bioplib/MathType.h"

#ifndef MATRIX 
#define MATRIX

int compute_covariance_matrix(double **x,double **cov,int numberOfPoints,int numberOfDimensions);

void find_mean(double **points,double *centroid,int numberOfPoints,int numberOfDimensions);

#endif
