/************************************************************************/
/**

   Program:    abpackingangle
   \file       abpackingangle.c
   
   \version    V1.3
   \date       03.10.16
   \brief      Calculate the VH/VL packing angle for an antibody Fv
   
   \copyright  (c) UCL / Abhi Raghavan and Dr. Andrew C. R. Martin 2007-16
   \author     Dr. Abhi Raghavan and Dr. Andrew C. R. Martin
   \par
               Institute of Structural & Molecular Biology,
               University College,
               Gower Street,
               London.
               WC1E 6BT.
   \par
               andrew@bioinf.org.uk
               andrew.martin@ucl.ac.uk
               
**************************************************************************

   This program is not in the public domain, but it may be copied
   according to the conditions laid out in the accompanying file
   COPYING.DOC

   The code may be modified as required, but any modifications must be
   documented so that the person responsible can be identified.

   The code may not be sold commercially or included as part of a 
   commercial product except as described in the file COPYING.DOC.

**************************************************************************

   Description:
   ============

**************************************************************************

   Usage:
   ======

**************************************************************************

   Revision History:
   =================
   V1.0   02.03.07 Original version By: Abhi
   V1.1   02.06.16 Split from Abhi's code base  By: ACRM
   V1.2   03.10.16 Some tidying up and changed to new Bioplib routines
                   builddist config
   V1.3   03.10.16 Full tidy up
   V1.4   20.10.16 Removed all global variables
                   Changed to my style command line parser

*************************************************************************/
/* Includes
*/
#include <math.h>
#include <string.h>
#include <strings.h>
#include <unistd.h>
#include <ctype.h>
#include <stdlib.h>

#include "bioplib/angle.h"
#include "bioplib/MathUtil.h"
#include "bioplib/SysDefs.h"
#include "bioplib/pdb.h"
#include "matrix.h"
#include "regression.h"

/************************************************************************/
/* Defines and macros
*/
#define MAXPOINTS    64
#define MAXFILENAME 160
#define MAXPDBCODE    8

/************************************************************************/
/* Globals
*/


/************************************************************************/
/* Prototypes
*/
int get_number_of_constant_residue_positions(char **constantPositions);
void free_pointers(void);
void create_constant_positions_list(
   char **lightChainConstantPositionsList,
   int *numberOfConstantLightChainPositions,
   char **heavyChainConstantPositionsList,
   int *numberOfConstantHeavyChainPositions);
int read_coordinates(char **chainConstantPositionsList,
                     REAL **coordArray, PDB *firstPointer);
REAL calculate_cosine_from_vectors(REAL *eigenVector1,  
                                     REAL *eigenVector2,  
                                     int numberOfDimensions,
                                     BOOL isUnitVector);
void draw_regression_line(REAL **coordinates,
                          REAL *eigenVector,
                          int numberOfPoints,
                          char *chainLabel,
                          FILE *wfp);
BOOL plot(FILE *fp1, FILE *fp2,
          FILE *wfp, BOOL displayOutputFlag, BOOL DisplayStatsFlag,
          char **lightChainConstantPositionsList,
          char **heavyChainConstantPositionsList);
void Usage(void);
BOOL ParseCommandLine(int argc, char **argv, char *lightFile, 
                      char *heavyFile, char *pdbCode,
                      char *vecFile, BOOL *debug);


/************************************************************************/
int main(int argc,char **argv)
{
   FILE *wfp=NULL, *fp1=NULL, *fp2=NULL;
   char LightChainFilename[MAXFILENAME],
        HeavyChainFilename[MAXFILENAME],
        outputFilename[MAXFILENAME],
        PDBCode[MAXPDBCODE],
        **lightChainConstantPositionsList=NULL,
        **heavyChainConstantPositionsList=NULL;

   BOOL DisplayOutputFlag = FALSE,
        DisplayStatsFlag  = FALSE;

   if(argc < 5)
   {
      Usage();
      return 0;
   }

   if(!ParseCommandLine(argc, argv, LightChainFilename, 
                        HeavyChainFilename, PDBCode, outputFilename,
                        &DisplayStatsFlag))
   {
      Usage();
      return 0;
   }

   if(outputFilename[0])
   {
      DisplayOutputFlag=TRUE;
   }
   

   if( access(LightChainFilename,R_OK) )
   {
      printf("\n Light chain file \"%s\" does not exist. Aborting program\n\n",
             LightChainFilename);
      return 1;
   }

   if( access(HeavyChainFilename,R_OK) )
   {
      printf("\n Heavy chain file \"%s\" does not exist. Aborting program\n\n",
             HeavyChainFilename);
      return 1;
   }

   lightChainConstantPositionsList = 
      (char **)malloc(MAXPOINTS * sizeof(char *));
   heavyChainConstantPositionsList = 
      (char **)malloc(MAXPOINTS * sizeof(char *));

   create_constant_positions_list(lightChainConstantPositionsList,NULL,
                                  heavyChainConstantPositionsList,NULL);

   if(PDBCode[0])
      printf("PDB Code: %s\n",PDBCode);

   if(DisplayOutputFlag)
   {
      if((wfp=fopen(outputFilename,"w"))==NULL)
      {
         fprintf(stderr,"Error (abpackingangle) - Can't write output file: %s\n",
                 outputFilename);
         return(1);
      }
      
   }

   if((fp1=fopen(LightChainFilename, "r"))==NULL)
   {
      fprintf(stderr,"Error (abpackingangle) - Can't read PDB file: %s\n",
              LightChainFilename);
      return(1);
   }
   
   if((fp2=fopen(HeavyChainFilename, "r"))==NULL)
   {
      fprintf(stderr,"Error (abpackingangle) - Can't read PDB file: %s\n",
              HeavyChainFilename);
      return(1);
   }
   

   if(!plot(fp1, fp2, wfp,
            DisplayOutputFlag, DisplayStatsFlag,
            lightChainConstantPositionsList,
            heavyChainConstantPositionsList))
   {
      return(1);
   }

   return 0;
}


/************************************************************************/
/*>int get_number_of_constant_residue_positions(char **constantPositions)
   ----------------------------------------------------------------------
*//**
   This function finds the number of constant residue positions
   defined in an array.  It returns the number of constant residue
   positions.  
*/
int get_number_of_constant_residue_positions(char **constantPositions)
{
   int i=0;

   while(constantPositions[i][0] != '0')
      i++;

   return i;

}


/************************************************************************/
void create_constant_positions_list(
   char **lightChainConstantPositionsList,
   int *numberOfConstantLightChainPositions,
   char **heavyChainConstantPositionsList,
   int *numberOfConstantHeavyChainPositions)
{
   int i=0;

   /* The constant residue positions are as follows:
      Light chain positions: 
         ["L35","L36","L37","L38"] and ["L85","L86","L87","L88"]
      Heavy chain positions: 
         ["H36","H37","H38","H39"] and ["H89","H90","H91","H92"]
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
   lightChainConstantPositionsList[8][0]='0';
   heavyChainConstantPositionsList[8][0]='0';

   /* Set the number of positions                                       */
   if(numberOfConstantLightChainPositions)
      *numberOfConstantLightChainPositions=8;

   if(numberOfConstantHeavyChainPositions)
      *numberOfConstantHeavyChainPositions=8;

} /* End of function "create_constant_positions_list". */


/************************************************************************/
/*>REAL calculate_cosine_from_vectors(REAL *vector1,
                                        REAL *vector2,
                                        int numberOfDimensions,
                                        BOOL isUnitVector)
   ------------------------------------------------------------
*//**
 \param[in] *vector1            First Vector
 \param[in] *vector2            Second Vector
 \param[in] numberOfDimensions  Number of dimensions
 \param[in] isUnitVector        Whether the Vectors are unit vectors.
*/
REAL calculate_cosine_from_vectors(REAL *vector1,
                                     REAL *vector2,
                                     int numberOfDimensions,
                                     BOOL isUnitVector)
{
   REAL numerator=0,
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

      return ( (REAL)numerator/denominator);
   }

   return numerator;

}


/************************************************************************/
int verify(REAL *projection, REAL *centroid, REAL *vector,
           char *string)
{
   REAL xc=0,
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


/************************************************************************/
REAL calculate_torsion_angle(REAL *lightChainVector,
                             REAL *lightChainCentroid,
                             REAL *lightChainPointToBeProjected,
                             REAL *heavyChainVector,
                             REAL *heavyChainCentroid,
                             REAL *heavyChainPointToBeProjected,
                             BOOL DisplayStatsFlag)
{
   /* Step 1: Declare variables to be used in the function. */

   REAL point[3],
          *lightChainProjection,
          *heavyChainProjection;

   int i=0;

   REAL torsionAngle=0;

   /* Step 2: Find two points on the light chain regression line. One
              of the points (default) is the centroid. The other point
              can be calculated the following way.
   */

   point[0]=lightChainCentroid[0] + (100 * lightChainVector[0]);
   point[1]=lightChainCentroid[1] + (100 * lightChainVector[1]);
   point[2]=lightChainCentroid[2] + (100 * lightChainVector[2]);

   /* Step 3: Find the point of projection on the light chain. We use
              the function blPointLineDistance to do this. The
              function syntax is as follows:

              REAL blPointLineDistance(REAL Px, REAL Py, REAL Pz,
                                       REAL P1x, REAL P1y, REAL P1z,
                                       REAL P2x, REAL P2y, REAL P2z,
                                       REAL *Rx, REAL *Ry, REAL *Rz,
                                       REAL *frac)

              (P1x,P1y,P1z) and (P2x,P2y,P2z) are two points on the
              line. Point (Px,Py,Pz) is to be projected onto this
              line.
   */

   lightChainProjection=(REAL *)malloc(3 * sizeof(REAL));

   blPointLineDistance(lightChainPointToBeProjected[0],
                       lightChainPointToBeProjected[1],
                       lightChainPointToBeProjected[2],
                       lightChainCentroid[0],
                       lightChainCentroid[1],
                       lightChainCentroid[2],
                       point[0],
                       point[1],
                       point[2],
                       &lightChainProjection[0],
                       &lightChainProjection[1],
                       &lightChainProjection[2],
                       NULL);

   verify(lightChainProjection, lightChainCentroid,
          lightChainVector,"Light");

   /* Step 4: Perform a similar procedure for the heavy chain. */
   point[0]=heavyChainCentroid[0] + (100 * heavyChainVector[0]);
   point[1]=heavyChainCentroid[1] + (100 * heavyChainVector[1]);
   point[2]=heavyChainCentroid[2] + (100 * heavyChainVector[2]);

   heavyChainProjection=(REAL *)malloc(3 * sizeof(REAL));

   blPointLineDistance(heavyChainPointToBeProjected[0],
                       heavyChainPointToBeProjected[1],
                       heavyChainPointToBeProjected[2],
                       heavyChainCentroid[0],
                       heavyChainCentroid[1],
                       heavyChainCentroid[2],
                       point[0],point[1],point[2],
                       &heavyChainProjection[0],
                       &heavyChainProjection[1],
                       &heavyChainProjection[2],
                       NULL);

   verify(heavyChainProjection, heavyChainCentroid,
          heavyChainVector, "Heavy");


   if(DisplayStatsFlag)
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

   /* Step 5: Now that the projected points have been found, find the
              torsion angle using the function "blPhi". The format of
              the function is as given below:

              REAL blPhi(REAL xi,
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

   torsionAngle=blPhi(lightChainProjection[0],
                      lightChainProjection[1],
                      lightChainProjection[2],
                      lightChainCentroid[0],
                      lightChainCentroid[1],
                      lightChainCentroid[2],
                      heavyChainCentroid[0],
                      heavyChainCentroid[1],
                      heavyChainCentroid[2],
                      heavyChainProjection[0],
                      heavyChainProjection[1],
                      heavyChainProjection[2]);

   /* Step 6: Return torsion angle to calling function. */

   return torsionAngle;

}




/************************************************************************/
/* void Usage(): Shows the usage of the program
*/

void Usage()
{
   printf("\n Usage: abpackingangle <Parameters>\n");
   printf("\n Parameters are:\n");
   printf("\n 1. -l <Light chain PDB file>");
   printf("\n 2. -h <Heavy chain PDB file>");
   printf("\n 3. -pdb <PDB Code (Optional parameter)>");
   printf("\n 4. -output <Name of Output file - PDB format file of light and heavy chain with best fit line (Optional parameter)>");
   printf("\n 4. -stats (Use this option if you would like to display important values calculated by the program)");
   printf("\n\n");

} /* End of function "Usage". */



/************************************************************************/
/*>BOOL ParseCommandLine(int argc, char **argv, char *lightFile, 
                         char *heavyFile, char *pdbCode, char *vecFile,
                         BOOL *debug)
   --------------------------------------------------------------------
*//**


*/
BOOL ParseCommandLine(int argc, char **argv, char *lightFile, 
                      char *heavyFile, char *pdbCode, char *vecFile,
                      BOOL *debug)
{
   lightFile[0] = '\0';
   heavyFile[0] = '\0';
   pdbCode[0]   = '\0';
   vecFile[0]   = '\0';
   

   argc--;
   argv++;
   
   while(argc)
   {
      if(argv[0][0] == '-')
      {
         switch(argv[0][1])
         {
/*
         case 'h':
            return(FALSE);
            break;
*/
         case 'l':
            argc--;
            argv++;
            if(!argc)
               return(FALSE);
            strncpy(lightFile, argv[0], MAXFILENAME);
            break;
         case 'h':
            argc--;
            argv++;
            if(!argc)
               return(FALSE);
            strncpy(heavyFile, argv[0], MAXFILENAME);
            break;
         case 'p':
            argc--;
            argv++;
            if(!argc)
               return(FALSE);
            strncpy(pdbCode, argv[0], MAXPDBCODE);
            break;
         case 'o':
            argc--;
            argv++;
            if(!argc)
               return(FALSE);
            strncpy(vecFile, argv[0], MAXFILENAME);
            break;
         case 's':
         case 'd':
            *debug = TRUE;
            break;
         default:
            return(FALSE);
         }
      }
      argc--;
      argv++;
   }
   
   return(TRUE);
}



/************************************************************************/
BOOL plot(FILE *fp1, FILE *fp2, FILE *wfp,
          BOOL displayOutputFlag, BOOL DisplayStatsFlag,
          char **lightChainConstantPositionsList,
          char **heavyChainConstantPositionsList)
{
   /* Step 1: Declare variables to be used in the function. */

   int lightChainNumberOfAtoms=0,
       nextPointIndex=-1,
       i=0,
       numberOfConstantLightChainPositions=0,
       numberOfConstantHeavyChainPositions=0;

   PDB *pdb=NULL;

   REAL **lightChainConstantPositionCoordinates=NULL,
      **heavyChainConstantPositionCoordinates=NULL,
      *lightEigenVector=NULL,
      *heavyEigenVector=NULL,
      **modifiedLightCoordinates=NULL,
      **modifiedHeavyCoordinates=NULL,
      *lightChainCentroid=NULL,
      *heavyChainCentroid=NULL,
      torsionAngle=0;

   /* Step 2: Read atoms from the light and heavy chain PDB files into
    * two linked lists */

   pdb=blReadPDB(fp1,&lightChainNumberOfAtoms);

   /* Step 3: Read the light and heavy chain constant position
              coordinates for CA atoms from the linked lists into
              arrays. Function:

      int read_coordinates(char **chainConstantPositionsList,
                           char **coordArray,PDB *firstPointer)
   */

   lightChainConstantPositionCoordinates = 
      (REAL **)malloc(MAXPOINTS * sizeof(REAL *));
   heavyChainConstantPositionCoordinates = 
      (REAL **)malloc(MAXPOINTS * sizeof(REAL *));

   if((numberOfConstantLightChainPositions = 
       read_coordinates(lightChainConstantPositionsList,
                        lightChainConstantPositionCoordinates,
                        pdb))==0)
   {
      return(FALSE);
   }
   

   if((numberOfConstantHeavyChainPositions = 
       read_coordinates(heavyChainConstantPositionsList,
                        heavyChainConstantPositionCoordinates,
                        pdb))==0)
   {
      return(FALSE);
   }

   /* Step 4: Find mid points of atoms that are structurally
              adjacent. This is done to find the best fit line.
   */

   /* Calculation of mid points for the Light chain. */

   modifiedLightCoordinates = 
      (REAL **)malloc( (numberOfConstantLightChainPositions/2) * 
                         sizeof(REAL *) );
   nextPointIndex=numberOfConstantLightChainPositions/2;

   for(i=0; i<numberOfConstantLightChainPositions/2; i++)
   {
      modifiedLightCoordinates[i]=(REAL *)malloc(3 * sizeof(REAL));

      modifiedLightCoordinates[i][0] = 
         (lightChainConstantPositionCoordinates[i][0] +
          lightChainConstantPositionCoordinates[i+nextPointIndex][0])/2;

      modifiedLightCoordinates[i][1] =
         (lightChainConstantPositionCoordinates[i][1] +
          lightChainConstantPositionCoordinates[i+nextPointIndex][1])/2;

      modifiedLightCoordinates[i][2] =
         (lightChainConstantPositionCoordinates[i][2] +
          lightChainConstantPositionCoordinates[i+nextPointIndex][2])/2;
   }

   /* Calculation of mid points for the heavy chain. */

   modifiedHeavyCoordinates = 
      (REAL **)malloc( (numberOfConstantHeavyChainPositions/2) * 
                         sizeof(REAL *) );
   nextPointIndex=numberOfConstantHeavyChainPositions/2;

   for(i=0; i<numberOfConstantHeavyChainPositions/2; i++)
   {
      modifiedHeavyCoordinates[i]=(REAL *)malloc(3 * sizeof(REAL));

      modifiedHeavyCoordinates[i][0] = 
         (heavyChainConstantPositionCoordinates[i][0] +
          heavyChainConstantPositionCoordinates[i+nextPointIndex][0])/2;
      
      modifiedHeavyCoordinates[i][1] =
         (heavyChainConstantPositionCoordinates[i][1] +
          heavyChainConstantPositionCoordinates[i+nextPointIndex][1])/2;

      modifiedHeavyCoordinates[i][2] =
         (heavyChainConstantPositionCoordinates[i][2] +
          heavyChainConstantPositionCoordinates[i+nextPointIndex][2])/2;
   }


   /* -------------- NEW BIT OF MODIFICATION ------------------ */

   /* Step 5: Call the routine to fit the least squares distant line
              for the points.  Format for function that performs this:

      void compute_best_fit_line(REAL **coordinates,
                                 int numberOfPoints,
                                 int numberOfDimensions,
                                 REAL *centroid,
                                 REAL *eigenVector)
   */

   lightEigenVector=(REAL *)malloc(3 * sizeof(REAL));
   heavyEigenVector=(REAL *)malloc(3 * sizeof(REAL));

   lightChainCentroid=(REAL *)malloc(3 * sizeof(REAL));
   heavyChainCentroid=(REAL *)malloc(3 * sizeof(REAL));
 
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

              Light chain: Light chain centroid and mid point of L35
              and L88.

              Heavy chain: Heavy chain centroid and mid point of L36
              and L92.

              We use the function "calculate_torsion_angle" whose
              syntax is as below:

              REAL calculate_torsion_angle(REAL *lightChainVector,
                                           REAL *lightChainCentroid,
                                           REAL *lightChainPointToBeProjected,
                                           REAL *heavyChainVector,
                                           REAL *heavyChainCentroid,
                                           REAL *heavyChainPointToBeProjected,
                                           BOOL   DisplayStatsFlag)
   */

   torsionAngle=calculate_torsion_angle(lightEigenVector,
                                        lightChainCentroid,
                                        modifiedLightCoordinates[0],
                                        heavyEigenVector,
                                        heavyChainCentroid,
                                        modifiedHeavyCoordinates[0], 
                                        DisplayStatsFlag);

   torsionAngle=(torsionAngle * 180)/PI;

   /* Step 7: Print the torsion angle and 3write the coordinates of
              the imaginary regression line into a PDB file (if
              required).
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

      blWritePDB(wfp,pdb);
   }

   printf("------------------------------------\n");

   /* Step 8: Release memory allocated to pointers and return to the
    * calling function 
    */

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

   FREELIST(pdb, PDB);

   return(TRUE);
}


/************************************************************************/
int read_coordinates(char **chainConstantPositionsList, REAL **coordArray,
                     PDB *pdb)
{
   int numberOfConstantChainPositions=0,
       i=0;

   PDB *p=NULL;

   numberOfConstantChainPositions =
      get_number_of_constant_residue_positions(chainConstantPositionsList);

   /* Get the X, Y, Z coordinates of the CA atoms in light and heavy
    * chain constant residue list.
    */

   for(i=0; i<numberOfConstantChainPositions; i++)
   {
      if((p = blFindResidueSpec(pdb, chainConstantPositionsList[i]))==NULL)
      {
         fprintf(stderr,"Error (abpackingangle): Residue not found (%s)\n",
                 chainConstantPositionsList[i]);
         return(0);
      }
      
      
      while(! strstr(p->atnam,"CA ") )
         p=p->next;

      if(! p)
         break;

      coordArray[i]=(REAL *)malloc(3 * sizeof(REAL));

      coordArray[i][0]=(REAL)p->x;
      coordArray[i][1]=(REAL)p->y;
      coordArray[i][2]=(REAL)p->z;
   }

   return i;

}


