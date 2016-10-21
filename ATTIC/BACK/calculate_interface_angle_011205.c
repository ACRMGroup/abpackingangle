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

FILE *fp=NULL;

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
}

int read_coordinates(char **chainConstantPositionsList,double **coordArray,PDB *firstPointer)
{
   int numberOfConstantChainPositions=0,
       i=0,
       atomNumber=0;

   char insertCode[8];

   PDB *p=NULL;

   /* Format of function:

      void get_number_of_constant_residue_positions(char **constantPositions)  
   */

   numberOfConstantChainPositions=get_number_of_constant_residue_positions(chainConstantPositionsList);

   for(i=0;i<numberOfConstantChainPositions;i++)
   {
      p=firstPointer;

      get_atom_number_and_insert_code(chainConstantPositionsList[i],&atomNumber,insertCode);

      while(p->resnum != atomNumber)
         p=p->next;

      while(! strstr(p->atnam,"CA ") )
         p=p->next;

      if(! p)
         break;

      coordArray[i]=(double *)malloc(3 * sizeof(double));

      coordArray[i][0]=(double)p->x;
      coordArray[i][1]=(double)p->y;
      coordArray[i][2]=(double)p->z;
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


double find_interface_angle(char *lightChainFilename,char *heavyChainFilename)
{
   FILE *fp=NULL;
   double angle=0;

   fp=fopen(lightChainFilename,"r");
   lightFirst=ReadPDB(fp,&lightChainNumberOfAtoms);
   fclose(fp);

   fp=fopen(heavyChainFilename,"r");
   heavyFirst=ReadPDB(fp,&heavyChainNumberOfAtoms);
   fclose(fp);

   /* First, read number of residues in the constant residues list for light and heavy chains */

   
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


   // Perform memory allocation for the light chain and heavy chain coordinates.

   lightCoordinates=(double **)malloc(MAXPOINTS * sizeof(double *));
   heavyCoordinates=(double **)malloc(MAXPOINTS * sizeof(double *));

   /* Read the Light and heavy chain constant position coordinates for CA atoms into an array. Format of function:

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


   /* Now, come to the part of choosing the eigen vector with the largest eigen value.

      We compare the eigen values for the light chain and store index of the largest
      Eigen value in the variable lightChainLargestEigenValueIndex.

      Similarly, we store the index of the largest Eigen value for the heavy chain
      in the variable heavyChainLargestEigenValueIndex.
   */

   /*
      double calculate_cosine(double **eigenVectors1,         // First Eigen vector
                              double *eigenValues1,           // First set of Eigen values
                              double **eigenVectors2,         // Second Eigen vector.
                              double *eigenValues2,           // Second set of Eigen values
                              int numberOfDimensions,         // Number of dimensions.
                              BOOL isLargestEigenValue,       // Flag to indicate type of eigen value for selection.
                              BOOL isUnitVector)              // Whether the eigen vectors are unit vectors.
   */

   angleCosine=calculate_cosine(lightEigenVectors,lightEigenValues,	// Light chain eigen vectors and values.
				heavyEigenVectors,heavyEigenValues,	// Heavy chain eigen vectors and values.
				3,					// Number of dimensions
				TRUE,TRUE);				// TRUE (Largest value), TRUE

   angle=(acos(angleCosine) * 180)/3.14159;

   return angle;

} // End of function "find_interface_angle".

int main(int argc,char **argv)
{
   /* Allocate memory for the light and heavy chain constant positions list.

      Light chain positions: "L35","L36","L37","L38","L85","L86","L87","L88"
      Heavy chain positions: "H36","H37","H38","H39","H89","H90","H91","H92"

      Copy these positions into the arrays.
   */

   if(argc < 5)
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

   angle=find_interface_angle(lightChainFilename,heavyChainFilename);

   if(displayPDBFlag)
      printf("PDB Code: %s\n",pdbCode);

   printf("Angle: %f\n--------------------------------------\n",angle);

   free_pointers();
}
