# include <stdlib.h>
# include "matrix.h"


/* void find_mean(REAL **points,REAL *centroid,int numberOfPoints,int numberOfDimensions):

   A function to find the centroid of a given set of points.

   points -> a 2D array containing the (x,y,z) coordinates of every point.
   centroid -> a 1D array containing the (x,y,z) coordinates of the centroid.
   numberOfPoints -> implicitly the number of rows in the array points.
   numberOfDimensions -> number of dimensions to the system - 2 if only x,y coordinates
			 are provided. 3 if x,y,z coordinates are provided.
*/

void find_mean(REAL **points,REAL *centroid,int numberOfPoints,int numberOfDimensions)
{
   int i=0,
       j=0;

   REAL total=0;

   for(j=0;j<numberOfDimensions;j++)
   {
      total=0;

      for(i=0;i<numberOfPoints;i++)
      {
         total+=points[i][j];
      }

      centroid[j]=total/numberOfPoints;
   }
}


/* int compute_covariance_matrix(REAL **x,               // The array with X, Y, Z coordinates.
                              REAL **cov,             // Covariance matrix that will be created.
                              int numberOfPoints,       // Number of points used to calculate matrix.
                              int numberOfDimensions)   // Number of dimensions (3 for points in 3D space).

   Function to find the covariance matrix for a given set of points
*/

int compute_covariance_matrix(REAL **x,		/* The array with X, Y, Z coordinates. */
			      REAL **cov,		/* Covariance matrix that will be created. */
			      int numberOfPoints,	/* Number of points used to calculate matrix. */
			      int numberOfDimensions)   /* Number of dimensions (3 for points in 3D space). */
{
   int i=0,
       j=0,
       start=0;

   REAL *centroid,
	  total=0;

   centroid=(REAL *)malloc(numberOfDimensions * sizeof(REAL));

   /* Step 1: Calculate the centroid of the given set of points. The centroid is stored as
	      the last element of the respective row. For example,

		X-component of the centroid -> x[numberOfPoints][0].
		Y-component of the centroid -> x[numberOfPoints][1].
		Z-component of the centroid -> x[numberOfPoints][2]
   */

   find_mean(x,centroid,numberOfPoints,numberOfDimensions);

   /* Step 2: Calculate the covariance matrix for the set of points. This is given by the
	      formula.

	      For i=0 to 3 (exclude 3)

		 For j=i to 3 (exclude 3)

		   Total=0

		   For start=0 to numberOfPoints (exclude numberOfPoints)

			Total+=( x[start][i] - centroid[i] ) *
			       ( x[start][j] - centroid[j] )

                               // End of start for loop 

		    Covariance(i,j)=Total/numberOfPoints.

                       // End of j for loop

                       // End of i for loop. 
   */

   for(i=0;i<numberOfDimensions;i++)
   {
      for(j=i;j<numberOfDimensions;j++)
      {
	 total=0;

	 for(start=0;start<numberOfPoints;start++)
	 {
	    total+=( (x[start][i] - centroid[i]) * (x[start][j] - centroid[j]) );
	 }

	 /*cov[i][j]=(REAL)total/(numberOfPoints); */
	 /*cov[j][i]=(REAL)total/(numberOfPoints); */

	 cov[i][j]=(REAL)total;
	 cov[j][i]=(REAL)total;
      }
   }

   free(centroid);

   return 1;
}
