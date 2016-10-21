# include <string.h>
# include <unistd.h>
# include "bioplib/pdb.h"
# include "tgmath.h"
# include "math.h"
# include "matrix.h"

# define MAXPOINTS 32

PDB *lightFirst=NULL,
    *heavyFirst=NULL;

char lightChainFilename[100],
     heavyChainFilename[100];

int lightChainNumberOfAtoms=0,
    heavyChainNumberOfAtoms=0;

char **lightChainConstantPositionsList=NULL,
     **heavyChainConstantPositionsList=NULL;

double **lightCoordinates=NULL,
       **heavyCoordinates=NULL;

double **lightEigenVectors=NULL,
       **heavyEigenVectors=NULL,
       *lightEigenValues=NULL,
       *heavyEigenValues=NULL;

double **lightCovarianceMatrix=NULL,
       **heavyCovarianceMatrix=NULL;

FILE *fp=NULL,
     *wfp=NULL;

int i=0,
    atomNumber=0,
    numberOfConstantLightChainPositions=0,
    numberOfConstantHeavyChainPositions=0;

BOOL displayPDBFlag=FALSE;

char pdbCode[8];

int lightChainLargestEigenValueIndex=0,
    heavyChainLargestEigenValueIndex=0;

double tempLightEigen=-1000,
       tempHeavyEigen=-1000;

double numerator=0,
       denominator=-1,
       den1=0,
       den2=0;

double angleCosine=0,
       angle=0;

char insertCode[8];


/* ----------------------- FUNCTION DEFINITION SECTION -------------------------- */


/* void get_atom_number_and_insert_code(char *constantPosition,int *atomNumber,char *insertCode):

   This function parses the constant position code into the constituent position in integer and
   insert code in string.

   For example, "L27A" is split into:

   *atomNumber=27
   *insertCode=A
*/

void get_atom_number_and_insert_code(char *constantPosition,int *atomNumber,char *insertCode)
{
   char *p=NULL,
	localConstantPosition[8],
	atomNumberString[8];

   int atomIndex=0,
       insertCodeIndex=0;

   strcpy(localConstantPosition,constantPosition);

   p=localConstantPosition;
   p++; // If the constant position is L35, skip the L.

   for(;*p != '\0';p++)
   {
      if( isdigit(*p) )
	 atomNumberString[atomIndex++]=*p;
      else
      if( isalpha(*p) )
	 insertCode[insertCodeIndex++]=*p;
   }

   atomNumberString[atomIndex]='\0';
   insertCode[insertCodeIndex]='\0';

   *atomNumber=atoi(atomNumberString);

} // End of function "get_atom_number_and_insert_code".


/* int get_number_of_constant_residue_positions(char **constantPositions):

   This function finds the number of constant residue positions defined in an array.
   It returns the number of constant residue positions.
*/

int get_number_of_constant_residue_positions(char **constantPositions)
{
   int i=0;

   while(constantPositions[i][0] != '0')
      i++;

   return i;

} // End of function "get_number_of_constant_residue_positions".


/* void free_pointers():

   This function releases memory allocated to pointers.
*/

void free_pointers()
{
   /* The following pointers need to be freed up.

      PDB *lightFirst=NULL,
          *heavyFirst=NULL;

      double **lightCoordinates=NULL,
             **heavyCoordinates=NULL;

      double **lightEigenVectors=NULL,
             **heavyEigenVectors=NULL,
             *lightEigenValues=NULL,
             *heavyEigenValues=NULL;

      double **lightCovarianceMatrix=NULL,
             **heavyCovarianceMatrix=NULL;
   */

   PDB *p=NULL,
       *prev=NULL;

   int i=0;

   // First, the PDB list pointers.

   p=lightFirst;

   while(p)
   {
      prev=p;
      p=p->next;
      free(prev);
   }

   p=heavyFirst;

   while(p)
   {
      prev=p;
      p=p->next;
      free(prev);
   }

   // Now the remaining pointers

   for(i=0;i<3;i++)
   {
      free(lightCoordinates[i]);
      free(heavyCoordinates[i]);

      free(lightEigenVectors[i]);
      free(heavyEigenVectors[i]);

      free(lightCovarianceMatrix[i]);
      free(heavyCovarianceMatrix[i]);

      free(lightChainConstantPositionsList[i]);
      free(heavyChainConstantPositionsList[i]);
   }

   for(i=3;i<numberOfConstantLightChainPositions;i++)
   {
      free(lightCoordinates[i]);
      free(lightChainConstantPositionsList[i]);
   }

   for(i=3;i<numberOfConstantHeavyChainPositions;i++)
   {
      free(heavyCoordinates[i]);
      free(heavyChainConstantPositionsList[i]);
   }

   free(lightChainConstantPositionsList);
   free(heavyChainConstantPositionsList);

   free(lightCoordinates);
   free(heavyCoordinates);

   free(lightEigenVectors);
   free(heavyEigenVectors);

   free(lightEigenValues);
   free(heavyEigenValues);

   free(lightCovarianceMatrix);
   free(heavyCovarianceMatrix);

} // End of function "free_pointers".


/* void Usage(): Shows the usage of the program
*/

void Usage()
{
   printf("\n Usage: ./calculate_interface_angle.exe <Parameters>\n");
   printf("\n Parameters are:\n");
   printf("\n 1. -l <Light chain PDB file>");
   printf("\n 2. -h <Heavy chain PDB file>");
   printf("\n 3. -pdb <PDB Code (Optional parameter)");
   printf("\n\n");

} // End of function "Usage".


/* BOOL parse_command_line_parameters(int numberOfParam,char **param):

   This function parses for command line parameters. TRUE is returned if all
   the parameters are successfully read into the program, else FALSE is returned.
*/


BOOL parse_command_line_parameters(int numberOfParam,char **param)
{
   int i=1;

   while(i < numberOfParam)
   {
      if(! strcasecmp(param[i],"-l") )
      {
	 strcpy(lightChainFilename,param[i+1]);
	 i+=2;
	 continue;
      }
      else
      if(! strcasecmp(param[i],"-h") )
      {
	 strcpy(heavyChainFilename,param[i+1]);
	 i+=2;
	 continue;
      }
      else
      if(! strcasecmp(param[i],"-pdb") )
      {
	 displayPDBFlag=TRUE;
	 strcpy(pdbCode,param[i+1]);
	 i+=2;
	 continue;
      }
      else
	 return FALSE;
   }

   return TRUE;

} // End of function "parse_command_line_parameters".


void create_constant_positions_list(char **lightChainConstantPositionsList,char **heavyChainConstantPositionsList)
{
   int i=0;

   /*

   for(i=0;i<9;i++)
   {
      lightChainConstantPositionsList[i]=(char *)malloc(8 * sizeof(char));
      heavyChainConstantPositionsList[i]=(char *)malloc(8 * sizeof(char));
   }

   for(i=35;i<=38;i++)
   {
      sprintf(lightChainConstantPositionsList[i-35],"L%d",i);
      sprintf(heavyChainConstantPositionsList[i-35],"H%d",(i+1));
   }

   for(i=85;i<=88;i++)
   {
      sprintf(lightChainConstantPositionsList[i-85+4],"L%d",i);
      sprintf(heavyChainConstantPositionsList[i-85+4],"H%d",(i+4));
   }

   */

   /* The constant residue positions are as follows:

      Light chain positions: "L35","L36","L37","L38","L85","L86","L87","L88"
      Heavy chain positions: "H36","H37","H38","H39","H89","H90","H91","H92"

      The positions L35 and L88 are structurally close. The pairs are as follows:

      L35 and L88
      L36 and L87
      L37 and L86
      L38 and L85

      H36 and H92
      H37 and H91
      H38 and H90
      H39 and H89

      We arrange the constant array positions in a way such that the i'th element and i+4'th
      elements are structurall close.
   */

   for(i=0;i<9;i++)
   {
      lightChainConstantPositionsList[i]=(char *)malloc(8 * sizeof(char));
      heavyChainConstantPositionsList[i]=(char *)malloc(8 * sizeof(char));
   }

   // Constant Light chain positions.

   strcpy(lightChainConstantPositionsList[0],"L35");
   strcpy(lightChainConstantPositionsList[1],"L36");
   strcpy(lightChainConstantPositionsList[2],"L37");
   strcpy(lightChainConstantPositionsList[3],"L38");

   strcpy(lightChainConstantPositionsList[4],"L88");
   strcpy(lightChainConstantPositionsList[5],"L87");
   strcpy(lightChainConstantPositionsList[6],"L86");
   strcpy(lightChainConstantPositionsList[7],"L85");

   // Constant Heavy chain positions.

   strcpy(heavyChainConstantPositionsList[0],"H36");
   strcpy(heavyChainConstantPositionsList[1],"H37");
   strcpy(heavyChainConstantPositionsList[2],"H38");
   strcpy(heavyChainConstantPositionsList[3],"H39");

   strcpy(heavyChainConstantPositionsList[4],"H92");
   strcpy(heavyChainConstantPositionsList[5],"H91");
   strcpy(heavyChainConstantPositionsList[6],"H90");
   strcpy(heavyChainConstantPositionsList[7],"H89");

   strcpy(lightChainConstantPositionsList[8],"0");
   strcpy(heavyChainConstantPositionsList[8],"0");
}

int read_coordinates(char **chainConstantPositionsList,double **coordArray,PDB *firstPointer)
{
   int numberOfConstantChainPositions=0,
       i=0,
       atomNumber1=0,
       atomNumber2=0;

   char insertCode1[8],
	insertCode2[8];

   PDB *p=NULL;

   /* Format of function:

      void get_number_of_constant_residue_positions(char **constantPositions)  
   */

   numberOfConstantChainPositions=get_number_of_constant_residue_positions(chainConstantPositionsList);
   numberOfConstantChainPositions-=4;

   /* Get the X, Y, Z coordinates of the CA atoms in light and heavy chain constant residue list.

      The format of the PDB structure is:

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

   for(i=0;i<numberOfConstantChainPositions;i++)
   {
      p=firstPointer;

      get_atom_number_and_insert_code(chainConstantPositionsList[i],&atomNumber1,insertCode1);
      get_atom_number_and_insert_code(chainConstantPositionsList[i],&atomNumber2,insertCode2);

      while(p->resnum != atomNumber1)
         p=p->next;

      while(! strstr(p->atnam,"CA ") )
         p=p->next;

      if(! p)
         break;

      coordArray[i]=(double *)malloc(3 * sizeof(double));

      coordArray[i][0]=(double)p->x;
      coordArray[i][1]=(double)p->y;
      coordArray[i][2]=(double)p->z;

      while(p->resnum != atomNumber2)
	 p=p->next;

      while(! strstr(p->atnam,"CA ") )
	 p=p->next;

      coordArray[
   }

   return numberOfConstantChainPositions;
}


void find_largest_smallest_indices_double(double *doubleArray,
				   	  int numberOfElements,
				   	  int *smallestValueIndex,
				   	  int *largestValueIndex)
{
   double smallestValue=1000000000,
	  largestValue=-1000000000;

   int i=0;

   for(i=0;i<numberOfElements;i++)
   {
      if(doubleArray[i] > largestValue)
      {
	 largestValue=doubleArray[i];

	 if(largestValueIndex)
	    *largestValueIndex=i;
      }

      if(doubleArray[i] < smallestValue)
      {
	 smallestValue=doubleArray[i];

	 if(smallestValueIndex)
	    *smallestValueIndex=i;
      }
   }

} // End of function "find_largest_smallest_indices_double".


double calculate_cosine(double **eigenVectors1,		// First Eigen vector
			double *eigenValues1,		// First set of Eigen values
			double **eigenVectors2,		// Second Eigen vector.
			double *eigenValues2,		// Second set of Eigen values
			int numberOfDimensions,		// Number of dimensions.
			BOOL isLargestEigenValue,	// Flag to indicate type of eigen value for selection.
			BOOL isUnitVector)		// Whether the eigen vectors are unit vectors.
{
   double numerator=0,
	  denominator=1,
	  den1=-1,
	  den2=-1;

   int i=0,
       index1=0,
       index2=0;

   /* Call the following function to find the index of the largest eigen values.

      void find_largest_smallest_indices_double(double *doubleArray,
                                                int numberOfElements,
                                                int *smallestValueIndex,
                                                int *largestValueIndex)
   */

   find_largest_smallest_indices_double(eigenValues1,numberOfDimensions,NULL,&index1);
   find_largest_smallest_indices_double(eigenValues2,numberOfDimensions,NULL,&index2);

   /* Find the Numerator of the cosine fraction */

   numerator=(eigenVectors1[index1][0] * eigenVectors2[index2][0])+
             (eigenVectors1[index1][1] * eigenVectors2[index2][1])+
             (eigenVectors1[index1][2] * eigenVectors2[index2][2]);

   if(! isUnitVector)
   {
      den1=(lightEigenVectors[index1][0]*lightEigenVectors[index1][0])+
           (lightEigenVectors[index1][1]*lightEigenVectors[index1][1])+
           (lightEigenVectors[index1][2]*lightEigenVectors[index1][2]);

      den2+=(heavyEigenVectors[index2][0]*heavyEigenVectors[index2][0])+
            (heavyEigenVectors[index2][1]*heavyEigenVectors[index2][1])+
            (heavyEigenVectors[index2][2]*heavyEigenVectors[index2][2]);

      denominator=den1 * den2;

      return (double)(numerator/denominator);
   }

   return numerator;
}


void plot_points(double **coordinates,
		 double *eigenVector,
		 int numberOfConstantChainPositions,
		 PDB *pdbListFirstPointer,
		 char *chainLabel,
		 FILE *wfp)
{
   PDB *p=NULL,
       *prev,
       *first=NULL;

   int smallestValueIndex=0,
       largestValueIndex=0;

   int kmin=0,
       kmax=0;

   int i=0;

   double xarray[10],
	  smallestX=0,
	  largestX=0;

   double *centroid;

   /* First find centroid of the light chain positions.

      void find_mean(double **points,double *centroid,int numberOfPoints,int numberOfDimensions)
   */

   centroid=(double *)malloc(3 * sizeof(double));

   find_mean(coordinates,centroid,numberOfConstantChainPositions,3);


   /* Use the function

      void find_largest_smallest_indices_double(double *doubleArray,
                                                int numberOfElements,
                                                int *smallestValueIndex,
                                                int *largestValueIndex)

      to get the indices of the smallest and largest X coordinates.
   */

   for(i=0;i<numberOfConstantChainPositions;i++)
   {
      xarray[i]=coordinates[i][0];
   }

   find_largest_smallest_indices_double(xarray,
					numberOfConstantChainPositions,
					&smallestValueIndex,
					&largestValueIndex);

   smallestX=xarray[smallestValueIndex];
   largestX=xarray[largestValueIndex];

   // We now ascertain the value of k so that points on the line of best fit may be plotted

   kmin=(int)((smallestX-centroid[0])/eigenVector[0]);
   kmax=(int)((largestX-centroid[0])/eigenVector[0]);

   if(kmin > kmax)
   {
      //kmin = -1 * kmin;
      //kmax = -1 * kmax;

      printf("\nKmin: %d\nKmax: %d",kmin,kmax);

      kmax=(-1 * kmax);
      kmin=0;

      printf("\nKmin: %d\nKmax: %d",kmin,kmax);
      getchar();
   }

   /* Now, create the new list of points using the values of k. PDB structure format is:

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

   for(i=kmin;i<kmax;i++)
   {
      p=(PDB *)malloc(1 * sizeof(PDB));

      strcpy(p->chain,chainLabel);

      p->x=centroid[0]+(i * eigenVector[0]);
      p->y=centroid[1]+(i * eigenVector[1]);
      p->z=centroid[2]+(i * eigenVector[2]);

      p->occ=1;
      p->bval=1;

      p->next=NULL;

      p->atnum=i-kmin+1;
      p->resnum=i-kmin+1;

      strcpy(p->record_type,"ATOM");
      strcpy(p->atnam,"O");
      strcpy(p->atnam_raw,"O");
      strcpy(p->resnam,"LYS");
      strcpy(p->insert,"");
      p->altpos=' ';

      if(prev)
	 prev->next=p;

      prev=p;

      if(first == NULL)
	 first = p;
   }

   // First, write the light chain.

   p=pdbListFirstPointer;

   while(p)
   {
      WritePDBRecord(wfp,p);
      p=p->next;
   }

   // Now write the PDB record of the vector.

   p=first;

   while(p)
   {
      WritePDBRecord(wfp,p);
      p=p->next;
   }
}



void plot(char *lightChainFilename,char *heavyChainFilename)
{
   FILE *fp=NULL;
   double angle=0;
   char outputFilename[100];

   int largestLightChainEigenValueIndex=0,
       largestHeavyChainEigenValueIndex=0;

   fp=fopen(lightChainFilename,"r");
   lightFirst=ReadPDB(fp,&lightChainNumberOfAtoms);
   fclose(fp);

   fp=fopen(heavyChainFilename,"r");
   heavyFirst=ReadPDB(fp,&heavyChainNumberOfAtoms);
   fclose(fp);

   // Perform memory allocation for the light chain and heavy chain coordinates.

   lightCoordinates=(double **)malloc(MAXPOINTS * sizeof(double *));
   heavyCoordinates=(double **)malloc(MAXPOINTS * sizeof(double *));

   /* Read the Light and heavy chain constant position coordinates for CA atoms into an array.

      Format of function:

      int read_coordinates(char **chainConstantPositionsList,char **coordArray,PDB *firstPointer)
   */

   numberOfConstantLightChainPositions=read_coordinates(lightChainConstantPositionsList,lightCoordinates,lightFirst);
   numberOfConstantHeavyChainPositions=read_coordinates(heavyChainConstantPositionsList,heavyCoordinates,heavyFirst);

   // Find number of constant light chain and heavy chain positions.

   /* Allocate space for the covariance matrices and compute them.

      The function that does the calculation of covariance matrices
      is declared as follows:

      int compute_covariance_matrix(double **x,               // The array with X, Y, Z coordinates.
                                    double **cov,             // Covariance matrix that will be created.
                                    int numberOfPoints,       // Number of points used to calculate matrix.
                                    int numberOfDimensions)   // Number of dimensions (3 for points in 3D space).
   */

   lightCovarianceMatrix=(double **)malloc(3 * sizeof(double *));
   heavyCovarianceMatrix=(double **)malloc(3 * sizeof(double *));

   for(i=0;i<3;i++)
   {
      lightCovarianceMatrix[i]=(double *)malloc(3 * sizeof(double));
      heavyCovarianceMatrix[i]=(double *)malloc(3 * sizeof(double));
   }

   compute_covariance_matrix(lightCoordinates,lightCovarianceMatrix,numberOfConstantLightChainPositions,3);
   compute_covariance_matrix(heavyCoordinates,heavyCovarianceMatrix,numberOfConstantHeavyChainPositions,3);

   /* Now, allocate space for the eigen values and vectors and calculate them by
      calling the appropriate function.

      The function that calculates Eigen values and vectors is declared as follows:

      int eigen(double **M,        // Input matrix
          	double **Vectors,  // eigen vector matrix -output
	        double *lambda,    // Output eigenvalues
        	int n)             // Input matrix dimension
   */

   lightEigenVectors=(double **)malloc(3 * sizeof(double *));
   heavyEigenVectors=(double **)malloc(3 * sizeof(double *));

   lightEigenValues=(double *)malloc(3 * sizeof(double));
   heavyEigenValues=(double *)malloc(3 * sizeof(double));

   for(i=0;i<3;i++)
   {
      lightEigenVectors[i]=(double *)malloc(3 * sizeof(double));
      heavyEigenVectors[i]=(double *)malloc(3 * sizeof(double));
   }

   eigen(lightCovarianceMatrix,lightEigenVectors,lightEigenValues,3);
   eigen(heavyCovarianceMatrix,heavyEigenVectors,heavyEigenValues,3);


   /* Now, come to the part of choosing the eigen vector with the largest eigen value.

      We compare the eigen values for the light chain and store index of the largest
      Eigen value in the variable lightChainLargestEigenValueIndex.

      Similarly, we store the index of the largest Eigen value for the heavy chain
      in the variable heavyChainLargestEigenValueIndex.
   */

   /*
      void find_largest_smallest_indices_double(double *doubleArray,
                                                int numberOfElements,
                                                int *smallestValueIndex,
                                                int *largestValueIndex)
   */

   find_largest_smallest_indices_double(lightEigenValues,3,NULL,&largestLightChainEigenValueIndex);
   find_largest_smallest_indices_double(heavyEigenValues,3,NULL,&largestHeavyChainEigenValueIndex);


   /* Call the function plot_points to write PDB records for the vector. Format of the function:

      void plot_points(double **coordinates,
                       double *eigenVector,
		       int numberOfConstantChainPositions,
                       PDB *pdbListFirstPointer,
                       FILE *wfp)
   */

   sprintf(outputFilename,"%s_.pdb",pdbCode);

   wfp=fopen(outputFilename,"w");

   plot_points(lightCoordinates,
	       lightEigenVectors[largestLightChainEigenValueIndex],
	       numberOfConstantLightChainPositions,
	       lightFirst,
	       "X",
	       wfp);

   /*
   plot_points(heavyCoordinates,
	       heavyEigenVectors[largestHeavyChainEigenValueIndex],
	       numberOfConstantHeavyChainPositions,
	       heavyFirst,
	       "Y",
	       wfp);
   */

   fclose(wfp);

} // End of function "find_interface_angle".


int main(int argc,char **argv)
{
   /* Allocate memory for the light and heavy chain constant positions list.

      Light chain positions: "L35","L36","L37","L38","L85","L86","L87","L88"
      Heavy chain positions: "H36","H37","H38","H39","H89","H90","H91","H92"

      Copy these positions into the arrays.
   */

   if(argc < 7)
   {
      Usage();
      exit(0);
   }

   if(! parse_command_line_parameters(argc,argv) )
   {
      Usage();
      exit(0);
   }

   if( access(lightChainFilename,R_OK) )
   {
      printf("\n Light chain file \"%s\" does not exist. Aborting program\n\n",lightChainFilename);
      exit(0);
   }

   if( access(heavyChainFilename,R_OK) )
   {
      printf("\n Heavy chain file \"%s\" does not exist. Aborting program\n\n",heavyChainFilename);
      exit(0);
   }

   lightChainConstantPositionsList=(char **)malloc(9 * sizeof(char *));
   heavyChainConstantPositionsList=(char **)malloc(9 * sizeof(char *));

   create_constant_positions_list(lightChainConstantPositionsList,heavyChainConstantPositionsList);

   strcpy(lightChainConstantPositionsList[8],"0");
   strcpy(heavyChainConstantPositionsList[8],"0");

   // Find the points along the least line vector and plot them.

   plot(lightChainFilename,heavyChainFilename);

   if(displayPDBFlag)
      printf("PDB Code: %s\n",pdbCode);

   /* Finally, free up all pointers.

      PDB *lightFirst=NULL,
          *heavyFirst=NULL;
      
      char lightChainFilename[100],
           heavyChainFilename[100];
      
      double **lightCoordinates=NULL,
             **heavyCoordinates=NULL;
      
      double **lightEigenVectors=NULL,
             **heavyEigenVectors=NULL,
             *lightEigenValues=NULL,
             *heavyEigenValues=NULL;
      
      double **lightCovarianceMatrix=NULL,
             **heavyCovarianceMatrix=NULL;
   */

   free_pointers();
}
