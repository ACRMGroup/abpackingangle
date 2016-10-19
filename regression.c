#include "bioplib/pdb.h"
#include "bioplib/eigen.h"
#include "arrays.h"
#include "matrix.h"

/* void DrawRegressionLine(REAL **coordinates,     // Coordinates of points used to calculate regression line
                             REAL *eigenVector,      // Eigen vector components of the regression line
                             int numberOfPoints,       // Number of points (number of rows in coordinates array).
                             char *chainLabel,         // Chain label when writing to PDB file.
                             FILE *wfp)                // File pointer to write into a PDB file.

   This function outputs an imaginary line in the form of PDB records. From the coordinates supplied,
   the centroid is calculated. The eigen vector components of the regression line are used to find
   other points on the line that passes through the centroid.

   The function returns no values.
*/


void DrawRegressionLine(REAL **coordinates,
                          REAL *eigenVector,
                          int numberOfPoints,
                          char *chainLabel,
                          FILE *wfp)
{
   /* ---------- Declare and define local variables ---------- */

   PDB *p=NULL,
       *prev,
       *first=NULL;

   int smallestValueIndex=0,
       largestValueIndex=0;

   int kmin=0,
       kmax=0;

   int i=0;

   REAL xarray[10],
          smallestX=0,
          largestX=0;

   REAL *centroid=NULL;


   /* Step 1: First find centroid of the light chain positions.

      void find_mean(REAL **points,REAL *centroid,int numberOfPoints,int numberOfDimensions)
   */

   centroid=(REAL *)malloc(3 * sizeof(REAL));

   find_mean(coordinates,centroid,numberOfPoints,3);


   /* Step 2: Use the function

      void find_largest_smallest_indices_REAL(REAL *REALArray,
                                                int numberOfElements,
                                                int *smallestValueIndex,
                                                int *largestValueIndex)

      to get the indices of the smallest and largest X coordinates.

      This will help find a range of values (k) over which points on the line
      can be plotted.
   */

   for(i=0;i<numberOfPoints;i++)
   {
      xarray[i]=coordinates[i][0];
   }

   find_largest_smallest_indices_REAL(xarray,
                                        numberOfPoints,
                                        &smallestValueIndex,
                                        &largestValueIndex);

   smallestX=xarray[smallestValueIndex];
   largestX=xarray[largestValueIndex];

   /* We now ascertain the value of k so that points on the line of best fit may be plotted */

   kmin=(int)((smallestX-centroid[0])/eigenVector[0]);
   kmax=(int)((largestX-centroid[0])/eigenVector[0]);

   if(kmin > kmax)
   {
      kmin = -1 * kmin;
      kmax = -1 * kmax;
   }

   /* Step 3: Create a new list for the points using the values of k.

      PDB structure format is:

      typedef struct pdb_entry
      {
         REAL x,y,z,occ,bval;
         struct pdb_entry *next;
         int  atnum;
         int  resnum;
         char record_type[8];
         char atnam[8];
         char atnam_raw[8];
         char resnam[8];
         char insert[8];
         char chain[8];
         char altpos;
      }  PDB;
   */

   p=NULL;
   prev=NULL;
   first=NULL;

   for(i=kmin;i<kmax;i++)
   {
      p=(PDB *)malloc(1 * sizeof(PDB));

      strcpy(p->chain,chainLabel);

      p->x=centroid[0]+(i * eigenVector[0]);
      p->y=centroid[1]+(i * eigenVector[1]);
      p->z=centroid[2]+(i * eigenVector[2]);

      p->occ=1;
      p->bval=1;

      if(prev)
         prev->next=p;

      p->next=NULL;

      p->atnum=i-kmin+1;
      p->resnum=i-kmin+1;

      p->atnum=i;
      p->resnum=i;

      strcpy(p->record_type,"ATOM");
      strcpy(p->atnam,"O");
      strcpy(p->atnam_raw,"O");
      strcpy(p->resnam,"LYS");
      strcpy(p->insert,"");
      p->altpos=' ';

      prev=p;

      if(first == NULL)
         first=p;
   }

   /* Step 4: Output PDB records for the regression line. Write to the file pointer wfp.
	      If wfp is NULL, then write to stdout.
	      Also, clear up space used for the list.
   */

   if(! wfp)
   {
      wfp=stdout;
   }

   p=first;

   while(p)
   {
      blWritePDBRecord(wfp,p);
      prev=p;
      p=p->next;
      free(prev);
   }

   first=NULL;

   /* Step 5: Release the space allocated to the centroid coordinates.
	      Exit from the function.
   */

   free(centroid);

   return;
}


void ComputeBestFitLine(REAL **coordinates,
                        int numberOfPoints,
                        int numberOfDimensions,
                        REAL *centroid,
                        REAL *eigenVector)
{
   /* Step 1: Declare the variables */

   REAL *eigenValues=NULL;

   REAL **covarianceMatrix=NULL,
	  **eigenVectorMatrix=NULL;

   int smallestEigenValueIndex=0,
       largestEigenValueIndex=0;

   int i=0;

   /* Step 2: Allocate memory for all the dynamic variables */

   eigenValues=(REAL *)malloc(numberOfDimensions * sizeof(REAL));

   covarianceMatrix=(REAL **)malloc(numberOfDimensions * sizeof(REAL *));
   eigenVectorMatrix=(REAL **)malloc(numberOfDimensions * sizeof(REAL *));

   for(i=0;i<numberOfDimensions;i++)
   {
      covarianceMatrix[i]=(REAL *)malloc(numberOfDimensions * sizeof(REAL));
      eigenVectorMatrix[i]=(REAL *)malloc(numberOfDimensions * sizeof(REAL));
   }

   /* Step 3: Find the centroid of the points. */

   find_mean(coordinates,centroid,numberOfPoints,numberOfDimensions);

   /* Step 4: Find the covariance matrix. */

   compute_covariance_matrix(coordinates,covarianceMatrix,numberOfPoints,numberOfDimensions);

   /* Step 5: Find the eigen values and vectors. */

   blEigen(covarianceMatrix,eigenVectorMatrix,eigenValues,numberOfDimensions);

   /* Step 6: Find the eigen vector that corresponds to the maximum eigen value. We use the function
	      "find_largest_smallest_indices_REAL" for this. The format for this function is:

	      void find_largest_smallest_indices_REAL(REAL *REALArray,
	 	                                        int numberOfElements,
                		                        int *smallestValueIndex,
                                		        int *largestValueIndex)
   */

   find_largest_smallest_indices_REAL(eigenValues,
					numberOfDimensions,
					&smallestEigenValueIndex,
					&largestEigenValueIndex);

   for(i=0;i<numberOfDimensions;i++)
   {
      eigenVector[i]=eigenVectorMatrix[i][largestEigenValueIndex];
   }

   /* Step 7: Final step, release all the memory allocated. */

   for(i=0;i<numberOfDimensions;i++)
   {
      free(covarianceMatrix[i]);
      free(eigenVectorMatrix[i]);
   }

   free(covarianceMatrix);
   free(eigenVectorMatrix);

   free(eigenValues);
}
