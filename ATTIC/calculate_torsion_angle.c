/* PREVIOUSLY, plot.c */

# include <string.h>
# include <strings.h>
# include <unistd.h>
# include <ctype.h>
# include <stdlib.h>
# include <math.h>

# include "bioplib/angle.h"
# include "bioplib/MathUtil.h"
# include "tgmath.h"
# include "matrix.h"
# include "regression.h"

# define MAXPOINTS 32

static char lightChainFilename[100],
            heavyChainFilename[100];

static BOOL displayPDBFlag=FALSE,
            displayOutputFlag=FALSE,
            displayStats=FALSE;

static char pdbCode[8];

static char **lightChainConstantPositionsList=NULL,
            **heavyChainConstantPositionsList=NULL;

static char outputFilename[100];


/* ---------------------------- LIST OF FUNCTIONS -------------------------------

void get_atom_number_and_insert_code(char *constantPosition,int *atomNumber,char *insertCode

int get_number_of_constant_residue_positions(char **constantPositions)

void free_pointers()

void create_constant_positions_list(char **lightChainConstantPositionsList,
                                    int *numberOfConstantLightChainPositions,
                                    char **heavyChainConstantPositionsList,
                                    int *numberOfConstantHeavyChainPositions)

int read_coordinates(char **chainConstantPositionsList,double **coordArray,PDB *firstPointer)

double calculate_cosine_from_vectors(double *eigenVector1,         // First Eigen vector
                                     double *eigenVector2,         // Second Eigen vector.
                                     int numberOfDimensions,       // Number of dimensions.
                                     BOOL isUnitVector)            // Whether the eigen vectors are unit vectors.

void draw_regression_line(double **coordinates,
			  double *eigenVector,
			  int numberOfPoints,
			  char *chainLabel,
			  FILE *wfp)

void plot(char *lightChainFilename,char *heavyChainFilename,FILE *wfp)

void Usage();

*/


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
   p++; /* Example - If position is L35, skip the L. */

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

} /* End of function "get_atom_number_and_insert_code". */


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

} /* End of function "get_number_of_constant_residue_positions". */



void create_constant_positions_list(char **lightChainConstantPositionsList,
				    int *numberOfConstantLightChainPositions,
				    char **heavyChainConstantPositionsList,
				    int *numberOfConstantHeavyChainPositions)
{
   int i=0;

   /* The constant residue positions are as follows:

      Light chain positions: ["L35","L36","L37","L38"] and ["L85","L86","L87","L88"]
      Heavy chain positions: ["H36","H37","H38","H39"] and ["H89","H90","H91","H92"]
   */

   for(i=0;i<9;i++)
   {
      lightChainConstantPositionsList[i]=(char *)malloc(8 * sizeof(char));
      heavyChainConstantPositionsList[i]=(char *)malloc(8 * sizeof(char));
   }

   /* Copy light chain constant positions. */

   strcpy(lightChainConstantPositionsList[0],"L35");
   strcpy(lightChainConstantPositionsList[1],"L36");
   strcpy(lightChainConstantPositionsList[2],"L37");
   strcpy(lightChainConstantPositionsList[3],"L38");

   strcpy(lightChainConstantPositionsList[4],"L88");
   strcpy(lightChainConstantPositionsList[5],"L87");
   strcpy(lightChainConstantPositionsList[6],"L86");
   strcpy(lightChainConstantPositionsList[7],"L85");

   /* Copy heavy chain constant positions. */

   strcpy(heavyChainConstantPositionsList[0],"H36");
   strcpy(heavyChainConstantPositionsList[1],"H37");
   strcpy(heavyChainConstantPositionsList[2],"H38");
   strcpy(heavyChainConstantPositionsList[3],"H39");

   strcpy(heavyChainConstantPositionsList[4],"H92");
   strcpy(heavyChainConstantPositionsList[5],"H91");
   strcpy(heavyChainConstantPositionsList[6],"H90");
   strcpy(heavyChainConstantPositionsList[7],"H89");

   /* Terminate the two lists. */

   if(numberOfConstantLightChainPositions)
      *numberOfConstantLightChainPositions=8;

   if(numberOfConstantHeavyChainPositions)
      *numberOfConstantHeavyChainPositions=8;

   lightChainConstantPositionsList[8][0]='0';
   heavyChainConstantPositionsList[8][0]='0';

} /* End of function "create_constant_positions_list". */


int read_coordinates(char **chainConstantPositionsList,double **coordArray,PDB *firstPointer)
{
   int numberOfConstantChainPositions=0,
       i=0,
       atomNumber1=0;

   char insertCode1[8];

   PDB *p=NULL;

   /* Format of function:

      void get_number_of_constant_residue_positions(char **constantPositions)  
   */

   numberOfConstantChainPositions=get_number_of_constant_residue_positions(chainConstantPositionsList);

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
   }

   return i;

} /* End of function "read_coordinates". */


double calculate_cosine_from_vectors(double *vector1,		/* First Vector */
		  		     double *vector2,		/* Second Vector. */
				     int numberOfDimensions,	/* Number of dimensions. */
				     BOOL isUnitVector)		/* Whether the Vectors are unit vectors. */
{
   double numerator=0,
	  denominator=1,
	  den1=-1,
	  den2=-1;

   int i=0;

   /* Find the Numerator of the cosine fraction */

   for(i=0;i<numberOfDimensions;i++)
   {
      numerator+=(vector1[i] * vector2[i]);
   }

   if(! isUnitVector)
   {
      den1=0;
      den2=0;

      for(i=0;i<numberOfDimensions;i++)
      {
	 den1+=(vector1[i] * vector1[i]);
	 den2+=(vector2[i] * vector2[i]);
      }

      den1 = sqrt(den1);
      den2 = sqrt(den2);

      denominator = den1 * den2;

      return ( (double)numerator/denominator);
   }

   return numerator;

} /* End of function "calculate_cosine_from_vectors". */


int verify(double *projection,double *centroid,double *vector,char *string)
{
   double xc=0,
          yc=0,
          zc=0,
          sq=0,
          delta=0,
	  den=0;

   xc=projection[0]-centroid[0];
   yc=projection[1]-centroid[1];
   zc=projection[2]-centroid[2];

   sq=(xc * xc) + (yc * yc) + (zc * zc);
   den=pow(sq,0.5);

   delta = (xc * xc)/(den * den) - (vector[0] * vector[0]);
   delta+= (yc * yc)/(den * den) - (vector[1] * vector[1]);
   delta+= (zc * zc)/(den * den) - (vector[2] * vector[2]);

   if(delta > 0.1)
      printf("\n Discrepancy in %s\n\n",string);

   return 1;
}


REAL calculate_torsion_angle(REAL *lightChainVector,
  	  		     REAL *lightChainCentroid,
  			     REAL *lightChainPointToBeProjected,
  			     REAL *heavyChainVector,
  			     REAL *heavyChainCentroid,
  			     REAL *heavyChainPointToBeProjected)
{
   /* Step 1: Declare variables to be used in the function. */

   REAL point[3],
	  *lightChainProjection,
	  *heavyChainProjection;

   int i=0;

   REAL torsionAngle=0;

   /* Step 2: Find two points on the light chain regression line. One of the points (default)
	      is the centroid. The other point can be calculated the following way.
   */

   point[0]=lightChainCentroid[0] + (100 * lightChainVector[0]);
   point[1]=lightChainCentroid[1] + (100 * lightChainVector[1]);
   point[2]=lightChainCentroid[2] + (100 * lightChainVector[2]);

   /* Step 3: Find the point of projection on the light chain. We use the function PointLineDistance
	      to do this. The function syntax is as follows:

	      REAL PointLineDistance(REAL Px, REAL Py, REAL Pz,
	                             REAL P1x, REAL P1y, REAL P1z,
	                             REAL P2x, REAL P2y, REAL P2z,
	                             REAL *Rx, REAL *Ry, REAL *Rz,
	                             REAL *frac)

	      (P1x,P1y,P1z) and (P2x,P2y,P2z) are two points on the line. Point (Px,Py,Pz) is
	      to be projected onto this line.
   */

   lightChainProjection=(REAL *)malloc(3 * sizeof(REAL));

   PointLineDistance(lightChainPointToBeProjected[0],lightChainPointToBeProjected[1],lightChainPointToBeProjected[2],
		     lightChainCentroid[0],lightChainCentroid[1],lightChainCentroid[2],
		     point[0],point[1],point[2],
		     &lightChainProjection[0],&lightChainProjection[1],&lightChainProjection[2],
		     NULL);

   verify(lightChainProjection,lightChainCentroid,lightChainVector,"Light");

   /* Step 4: Perform a similar procedure for the heavy chain. */

   point[0]=heavyChainCentroid[0] + (100 * heavyChainVector[0]);
   point[1]=heavyChainCentroid[1] + (100 * heavyChainVector[1]);
   point[2]=heavyChainCentroid[2] + (100 * heavyChainVector[2]);

   heavyChainProjection=(REAL *)malloc(3 * sizeof(REAL));

   PointLineDistance(heavyChainPointToBeProjected[0],heavyChainPointToBeProjected[1],heavyChainPointToBeProjected[2],
                     heavyChainCentroid[0],heavyChainCentroid[1],heavyChainCentroid[2],
                     point[0],point[1],point[2],
                     &heavyChainProjection[0],&heavyChainProjection[1],&heavyChainProjection[2],
                     NULL);

   verify(heavyChainProjection,heavyChainCentroid,heavyChainVector,"Heavy");


   if(displayStats)
   {
      printf("\nLight chain centroid:\t");

      for(i=0;i<3;i++)
      {
	 printf("\t%f",lightChainCentroid[i]);
      }

      printf("\nLight chain projection:\t");

      for(i=0;i<3;i++)
      {
	 printf("\t%f",lightChainProjection[i]);
      }

      printf("\nHeavy chain centroid:\t");

      for(i=0;i<3;i++)
      {
	 printf("\t%f",heavyChainCentroid[i]);
      }

      printf("\nHeavy chain projection:\t");

      for(i=0;i<3;i++)
      {
	 printf("\t%f",heavyChainProjection[i]);
      }

      printf("\n\n");
   }

   /* Step 5: Now that the projected points have been found, find the torsion angle using the
	      function "phi". The format of the function is as given below:

	      REAL phi(REAL xi,
	               REAL yi,
	               REAL zi,
	               REAL xj,
	               REAL yj,
	               REAL zj,
	               REAL xk,
	               REAL yk,
	               REAL zk,
	               REAL xl,
	               REAL yl,
		       REAL zl)
   */

   torsionAngle=phi(lightChainProjection[0],lightChainProjection[1],lightChainProjection[2],
		    lightChainCentroid[0],lightChainCentroid[1],lightChainCentroid[2],
		    heavyChainCentroid[0],heavyChainCentroid[1],heavyChainCentroid[2],
   		    heavyChainProjection[0],heavyChainProjection[1],heavyChainProjection[2]);

   /* Step 6: Return torsion angle to calling function. */

   return torsionAngle;

} /* end of function "calculate_torsion_angle". */



void plot(char *lightChainFilename,char *heavyChainFilename,FILE *wfp)
{
   /* Step 1: Declare variables to be used in the function. */

   FILE *fp=NULL;

   double **modifiedLightCoordinates=NULL,
	  **modifiedHeavyCoordinates=NULL;

   int lightChainNumberOfAtoms=0,
       heavyChainNumberOfAtoms=0;

   int nextPointIndex=-1;

   int i=0;

   PDB *lightFirst=NULL,
       *heavyFirst=NULL;

   PDB *p=NULL,
       *prev=NULL;

   double **lightChainConstantPositionCoordinates=NULL,
	  **heavyChainConstantPositionCoordinates=NULL;

   double *lightEigenVector=NULL,
	  *heavyEigenVector=NULL;

   int numberOfConstantLightChainPositions=0,
       numberOfConstantHeavyChainPositions=0;

   double *lightChainCentroid=NULL,
	  *heavyChainCentroid=NULL;

   double torsionAngle=0;

   /* Step 2: Read atoms from the light and heavy chain PDB files into two linked lists */

   fp=fopen(lightChainFilename,"r");
   lightFirst=ReadPDB(fp,&lightChainNumberOfAtoms);
   fclose(fp);

   fp=fopen(heavyChainFilename,"r");
   heavyFirst=ReadPDB(fp,&heavyChainNumberOfAtoms);
   fclose(fp);

   /* Step 3: Read the light and heavy chain constant position coordinates for CA atoms
	      from the linked lists into arrays. Function:

      int read_coordinates(char **chainConstantPositionsList,char **coordArray,PDB *firstPointer)
   */

   lightChainConstantPositionCoordinates=(double **)malloc(MAXPOINTS * sizeof(double *));
   heavyChainConstantPositionCoordinates=(double **)malloc(MAXPOINTS * sizeof(double *));

   numberOfConstantLightChainPositions=read_coordinates(lightChainConstantPositionsList,
							lightChainConstantPositionCoordinates,
							lightFirst);

   numberOfConstantHeavyChainPositions=read_coordinates(heavyChainConstantPositionsList,
							heavyChainConstantPositionCoordinates,
							heavyFirst);

   /* Step 4: Find mid points of atoms that are structurally adjacent. This is done to find
	      the best fit line.
   */

   /* Calculation of mid points for the Light chain. */

   modifiedLightCoordinates=(double  **)malloc( (numberOfConstantLightChainPositions/2) * sizeof(double *) );
   nextPointIndex=numberOfConstantLightChainPositions/2;

   for(i=0;i<numberOfConstantLightChainPositions/2;i++)
   {
      modifiedLightCoordinates[i]=(double *)malloc(3 * sizeof(double));

      modifiedLightCoordinates[i][0] = (lightChainConstantPositionCoordinates[i][0] +
					lightChainConstantPositionCoordinates[i+nextPointIndex][0])/2;

      modifiedLightCoordinates[i][1] = (lightChainConstantPositionCoordinates[i][1] +
					lightChainConstantPositionCoordinates[i+nextPointIndex][1])/2;

      modifiedLightCoordinates[i][2] = (lightChainConstantPositionCoordinates[i][2] +
					lightChainConstantPositionCoordinates[i+nextPointIndex][2])/2;
   }

   /* Calculation of mid points for the heavy chain. */

   modifiedHeavyCoordinates=(double **)malloc( (numberOfConstantHeavyChainPositions/2) * sizeof(double *) );
   nextPointIndex=numberOfConstantHeavyChainPositions/2;

   for(i=0;i<numberOfConstantHeavyChainPositions/2;i++)
   {
      modifiedHeavyCoordinates[i]=(double *)malloc(3 * sizeof(double));

      modifiedHeavyCoordinates[i][0] = (heavyChainConstantPositionCoordinates[i][0] +
					heavyChainConstantPositionCoordinates[i+nextPointIndex][0])/2;

      modifiedHeavyCoordinates[i][1] = (heavyChainConstantPositionCoordinates[i][1] +
					heavyChainConstantPositionCoordinates[i+nextPointIndex][1])/2;

      modifiedHeavyCoordinates[i][2] = (heavyChainConstantPositionCoordinates[i][2] +
					heavyChainConstantPositionCoordinates[i+nextPointIndex][2])/2;
   }


   /* -------------- NEW BIT OF MODIFICATION ------------------ */

   /* Step 5: Call the routine to fit the least squares distant line for the points.
              Format for function that performs this:

      void compute_best_fit_line(double **coordinates,
      				 int numberOfPoints,
      				 int numberOfDimensions,
      				 double *centroid,
      				 double *eigenVector)
   */

   lightEigenVector=(double *)malloc(3 * sizeof(double));
   heavyEigenVector=(double *)malloc(3 * sizeof(double));

   lightChainCentroid=(double *)malloc(3 * sizeof(double));
   heavyChainCentroid=(double *)malloc(3 * sizeof(double));
 
   compute_best_fit_line(modifiedLightCoordinates,
			 numberOfConstantLightChainPositions/2,
			 3,
			 lightChainCentroid,
			 lightEigenVector);

   compute_best_fit_line(modifiedHeavyCoordinates,
			 numberOfConstantHeavyChainPositions/2,
			 3,
			 heavyChainCentroid,
			 heavyEigenVector);

   /* Step 6: Calculate the torsion angle using the following points:

	      Light chain: Light chain centroid and mid point of L35 and L88.
	      Heavy chain: Heavy chain centroid and mid point of L36 and L92.

	      We use the function "calculate_torsion_angle" whose syntax is as below:

	      double calculate_torsion_angle(double *lightChainVector,
	                                     double *lightChainCentroid,
	                                     double *lightChainPointToBeProjected,
	                                     double *heavyChainVector,
	                                     double *heavyChainCentroid,
	                                     double *heavyChainPointToBeProjected)
   */

   torsionAngle=calculate_torsion_angle(lightEigenVector,
			                lightChainCentroid,
			                modifiedLightCoordinates[0],
			                heavyEigenVector,
			                heavyChainCentroid,
			                modifiedHeavyCoordinates[0]);

   torsionAngle=(torsionAngle * 180)/PI;

   /* Step 7: Print the torsion angle and 3write the coordinates of the imaginary regression line
	      into a PDB file (if required).
   */

   printf("Torsion angle: %f\n",torsionAngle);

   if(displayOutputFlag)
   {
      draw_regression_line(lightChainConstantPositionCoordinates,
                           lightEigenVector,
                           numberOfConstantLightChainPositions,
                           "X",
                           wfp);
   
      draw_regression_line(heavyChainConstantPositionCoordinates,
                           heavyEigenVector,
                           numberOfConstantHeavyChainPositions,
                           "Y",
                           wfp);

      WritePDB(wfp,lightFirst);
      WritePDB(wfp,heavyFirst);
   }

   printf("------------------------------------\n");

   /* Step 8: Release memory allocated to pointers and return to the calling function */

   for(i=0;i<numberOfConstantLightChainPositions/2;i++)
   {
      free(modifiedLightCoordinates[i]);
   }

   for(i=0;i<numberOfConstantHeavyChainPositions/2;i++)
   {
      free(modifiedHeavyCoordinates[i]);
   }

   for(i=0;i<numberOfConstantLightChainPositions;i++)
   {
      free(lightChainConstantPositionCoordinates[i]);
   }

   for(i=0;i<numberOfConstantHeavyChainPositions;i++)
   {
      free(heavyChainConstantPositionCoordinates[i]);
   }

   free(modifiedLightCoordinates);
   free(modifiedHeavyCoordinates);
   free(lightChainConstantPositionCoordinates);
   free(heavyChainConstantPositionCoordinates);

   free(lightEigenVector);
   free(heavyEigenVector);

   free(lightChainCentroid);
   free(heavyChainCentroid);

   p=lightFirst;
   prev=NULL;

   while(p)
   {
      prev=p;
      p=p->next;
      free(prev);
   }

   p=heavyFirst;
   prev=NULL;

   while(p)
   {
      prev=p;
      p=p->next;
      free(prev);
   }

   /* Step : Close the file pointer and exit from the function */

} /* End of function "plot". */


/* void Usage(): Shows the usage of the program
*/

void Usage()
{
   printf("\n Usage: ./calculate_interface_angle.exe <Parameters>\n");
   printf("\n Parameters are:\n");
   printf("\n 1. -l <Light chain PDB file>");
   printf("\n 2. -h <Heavy chain PDB file>");
   printf("\n 3. -pdb <PDB Code (Optional parameter)>");
   printf("\n 4. -output <Name of Output file - PDB format file of light and heavy chain with best fit line (Optional parameter)>");
   printf("\n 4. -stats (Use this option if you would like to display important values calculated by the program)");
   printf("\n\n");

} /* End of function "Usage". */


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
      if(! strcasecmp(param[i],"-output") )
      {
	 displayOutputFlag=TRUE;
	 strcpy(outputFilename,param[i+1]);
	 i+=2;
	 continue;
      }
      else
      if(! strcasecmp(param[i],"-stats") )
      {
	 displayStats=TRUE;
	 i+=1;
	 continue;
      }
      else
         return FALSE;
   }

   return TRUE;

} /* End of function "parse_command_line_parameters". */


int main(int argc,char **argv)
{
   /* Declare variables local to main. */

   FILE *wfp=NULL;

   if(argc < 5)
   {
      Usage();
      return 0;
   }

   if(! parse_command_line_parameters(argc,argv) )
   {
      Usage();
      return 0;
   }

   if( access(lightChainFilename,R_OK) )
   {
      printf("\n Light chain file \"%s\" does not exist. Aborting program\n\n",lightChainFilename);
      return 0;
   }

   if( access(heavyChainFilename,R_OK) )
   {
      printf("\n Heavy chain file \"%s\" does not exist. Aborting program\n\n",heavyChainFilename);
      return 0;
   }

   lightChainConstantPositionsList=(char **)malloc(MAXPOINTS * sizeof(char *));
   heavyChainConstantPositionsList=(char **)malloc(MAXPOINTS * sizeof(char *));

   create_constant_positions_list(lightChainConstantPositionsList,NULL,
				  heavyChainConstantPositionsList,NULL);

   if(displayPDBFlag)
      printf("PDB Code: %s\n",pdbCode);

   if(displayOutputFlag)
   {
      if(displayPDBFlag)
      {
	 wfp=fopen(outputFilename,"w");
      }
      else
      {
	 wfp=stdout;
      }
   }

   plot(lightChainFilename,heavyChainFilename,wfp);

   return 1;
}
