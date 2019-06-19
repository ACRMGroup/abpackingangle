/************************************************************************/
/**

   Program:    abpackingangle
   \file       abpackingangle.c
   
   \version    V2.1
   \date       19.06.19
   \brief      Calculate the VH/VL packing angle for an antibody Fv
   
   \copyright  (c) UCL / Abhi Raghavan and Dr. Andrew C. R. Martin 2007-19
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
   V2.0   19.10.16 Removed all global variables
                   Changed to my style command line parser
                   Now takes just one PDB file with both chains
                   PDB file now taken as parameter instead of switch
                   Added -q 
                   Simplified initialization of residue lists
                   Corrected usage message
                   Check all memory alocations
                   Use Bioplib routines where possible
                   More general cleanup
   V2.1   19.06.19 Added distance check between centroids and -d and -f
                   flags

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
#include "bioplib/macros.h"
#include "bioplib/MathUtil.h"
#include "bioplib/SysDefs.h"
#include "bioplib/pdb.h"
#include "bioplib/array.h"
#include "matrix.h"
#include "regression.h"

/************************************************************************/
/* Defines and macros
*/
#define MAXFILENAME 160
#define MAXPDBCODE  160
#define NUMRESPOS     8 /* Number of residue positions used for the 
                           vectors                                      */
#define MAXDIST      25 /* Highest distance seen with an OK (but odd)
                           structure is 19.636 (1MCO_1). 
                           Lowest distance seen with a clear wrong pairing
                           is 32.145 (3J42_2).                          */

/************************************************************************/
/* Globals
*/


/************************************************************************/
/* Prototypes
*/
int ReadCoordinates(char **keyResidueList,
                    REAL **coordArray,
                    PDB *firstPointer);
REAL CalculateCosineFromVectors(REAL *eigenVector1,  
                                REAL *eigenVector2,  
                                int  numberOfDimensions,
                                BOOL isUnitVector);
BOOL Plot(FILE *fpIn, FILE *fpOut,
          FILE *fpVecOut, BOOL displayOutputFlag, BOOL DisplayStatsFlag,
          BOOL quiet,
          char **lightKeyResidues,
          char **heavyKeyResidues,
          BOOL printDistance, BOOL force);
void Usage(void);
BOOL ParseCmdLine(int argc, char **argv, char *inFile, 
                  char *outFile, char *pdbCode,
                  char *vecFile, BOOL *verbose, BOOL *quiet,
                  BOOL *printDistance, BOOL *force);
void VerifyVector(REAL *projection, REAL *centroid, REAL *vector,
                  char *string);

/************************************************************************/
/*>int main(int argc, char **argv)
   -------------------------------
*//**
   Main program

-  19.10.16 Updated for V2.0  By: ACRM
*/
int main(int argc,char **argv)
{
   FILE *fpVecOut = NULL, 
        *fpIn     = stdin,
        *fpOut    = stdout;
   char inFilename[MAXFILENAME],
        outFilename[MAXFILENAME],
        vecFilename[MAXFILENAME],
        PDBCode[MAXPDBCODE];
   BOOL DisplayOutputFlag = FALSE,
        verbose           = FALSE,
        quiet             = FALSE,
        printDistance     = FALSE,
        force             = FALSE;

   char *lightKeyResidues[] = 
        {"L35", "L36", "L37", "L38", "L88", "L87", "L86", "L85", "\0"};
   char *heavyKeyResidues[] = 
        {"H36", "H37", "H38", "H39", "H92", "H91", "H90", "H89", "\0"};

   if(!ParseCmdLine(argc, argv, inFilename, outFilename, PDBCode, 
                    vecFilename, &verbose, &quiet, &printDistance, &force))
   {
      Usage();
      return 0;
   }

   if(vecFilename[0])
   {
      DisplayOutputFlag=TRUE;
   }
   
   /* Open files                                                        */
   if(!blOpenStdFiles(inFilename, outFilename, &fpIn, &fpOut))
   {
      fprintf(stderr,"Error (abpackingangle): Unable to open input or \
output file.\n");
      return(1);
   }
   if(DisplayOutputFlag)
   {
      if((fpVecOut=fopen(vecFilename,"w"))==NULL)
      {
         fprintf(stderr,"Error (abpackingangle): Can't write output \
vector PDB file: %s\n", vecFilename);
         return(1);
      }
   }

   if(PDBCode[0])
      fprintf(fpOut, "%s: ",PDBCode);

   if(!Plot(fpIn, fpOut, fpVecOut, DisplayOutputFlag, verbose, quiet,
            lightKeyResidues, heavyKeyResidues, printDistance, force))
   {
      return(1);
   }

   return 0;
}


/************************************************************************/
/*>REAL CalculateCosineFromVectors(REAL *vector1,
                                   REAL *vector2,
                                   int  numberOfDimensions,
                                   BOOL isUnitVector)
   ------------------------------------------------------------
*//**
 \param[in] *vector1            First Vector
 \param[in] *vector2            Second Vector
 \param[in] numberOfDimensions  Number of dimensions
 \param[in] isUnitVector        Whether the Vectors are unit vectors.
*/
REAL CalculateCosineFromVectors(REAL *vector1,
                                REAL *vector2,
                                int  numberOfDimensions,
                                BOOL isUnitVector)
{
   REAL numerator   =  0,
        denominator =  1,
        den1        = -1,
        den2        = -1;
   int  i = 0;

   /* Find the Numerator of the cosine fraction                         */
   for(i=0; i<numberOfDimensions; i++)
   {
      numerator += (vector1[i] * vector2[i]);
   }

   if(!isUnitVector)
   {
      den1 = 0;
      den2 = 0;

      for(i=0; i<numberOfDimensions; i++)
      {
         den1 += (vector1[i] * vector1[i]);
         den2 += (vector2[i] * vector2[i]);
      }

      den1 = sqrt(den1);
      den2 = sqrt(den2);

      denominator = den1 * den2;

      return (numerator/denominator);
   }

   return numerator;
}


/************************************************************************/
/*>void VerifyVector(REAL *projection, REAL *centroid, REAL *vector,
                     char *string)
   -----------------------------------------------------------------
*//**
*/
void VerifyVector(REAL *projection, REAL *centroid, REAL *vector, 
                  char *string)
{
   REAL xc    = 0,
        yc    = 0,
        zc    = 0,
        sq    = 0,
        delta = 0,
        den   = 0;

   xc  = projection[0]-centroid[0];
   yc  = projection[1]-centroid[1];
   zc  = projection[2]-centroid[2];
   sq  = (xc * xc) + (yc * yc) + (zc * zc);
   den = pow(sq,0.5);

   delta  = (xc * xc)/(den * den) - (vector[0] * vector[0]);
   delta += (yc * yc)/(den * den) - (vector[1] * vector[1]);
   delta += (zc * zc)/(den * den) - (vector[2] * vector[2]);

   if(delta > 0.1)
   {
      fprintf(stderr,"Warning (abpackingangle): Discrepancy in %s \
vector\n",string);
   }
}


/************************************************************************/
/*>REAL calculate_torsion_angle(REAL *lightVector,
                                REAL *lightCentroid,
                                REAL *lightPointToBeProjected,
                                REAL *heavyVector,
                                REAL *heavyCentroid,
                                REAL *heavyPointToBeProjected,
                                BOOL verbose)
   ------------------------------------------------------------
*//**
*/
REAL calculate_torsion_angle(REAL *lightVector,
                             REAL *lightCentroid,
                             REAL *lightPointToBeProjected,
                             REAL *heavyVector,
                             REAL *heavyCentroid,
                             REAL *heavyPointToBeProjected,
                             BOOL verbose)
{
   REAL point[3],
        lightProjection[3],
        heavyProjection[3],
        torsionAngle = 0;
   int  i            = 0;


   /* Find two points on the light chain regression line. One of the
      points (default) is the centroid. The other point can be
      calculated the following way.
   */
   point[0]=lightCentroid[0] + (100 * lightVector[0]);
   point[1]=lightCentroid[1] + (100 * lightVector[1]);
   point[2]=lightCentroid[2] + (100 * lightVector[2]);

   /* Find the point of projection on the light chain. We use the
      function blPointLineDistance to do this. The function syntax is
      as follows:

      REAL blPointLineDistance(REAL Px, REAL Py, REAL Pz,
                               REAL P1x, REAL P1y, REAL P1z,
                               REAL P2x, REAL P2y, REAL P2z,
                               REAL *Rx, REAL *Ry, REAL *Rz,
                               REAL *frac)

      (P1x,P1y,P1z) and (P2x,P2y,P2z) are two points on the
      line. Point (Px,Py,Pz) is to be projected onto this line.
   */
   blPointLineDistance(lightPointToBeProjected[0],
                       lightPointToBeProjected[1],
                       lightPointToBeProjected[2],
                       lightCentroid[0],
                       lightCentroid[1],
                       lightCentroid[2],
                       point[0],
                       point[1],
                       point[2],
                       &lightProjection[0],
                       &lightProjection[1],
                       &lightProjection[2],
                       NULL);

   VerifyVector(lightProjection, lightCentroid,
          lightVector,"Light");

   /* Do the same for the heavy chain.                                  */
   point[0] = heavyCentroid[0] + (100 * heavyVector[0]);
   point[1] = heavyCentroid[1] + (100 * heavyVector[1]);
   point[2] = heavyCentroid[2] + (100 * heavyVector[2]);

   blPointLineDistance(heavyPointToBeProjected[0],
                       heavyPointToBeProjected[1],
                       heavyPointToBeProjected[2],
                       heavyCentroid[0],
                       heavyCentroid[1],
                       heavyCentroid[2],
                       point[0],
                       point[1],
                       point[2],
                       &heavyProjection[0],
                       &heavyProjection[1],
                       &heavyProjection[2],
                       NULL);

   VerifyVector(heavyProjection, heavyCentroid,
          heavyVector, "Heavy");

   /* Print information if required                                     */
   if(verbose)
   {
      fprintf(stderr,"\nLight chain centroid:\t");

      for(i=0;i<3;i++)
         fprintf(stderr,"\t%f",lightCentroid[i]);

      fprintf(stderr,"\nLight chain projection:\t");

      for(i=0;i<3;i++)
         fprintf(stderr,"\t%f",lightProjection[i]);

      fprintf(stderr,"\nHeavy chain centroid:\t");

      for(i=0;i<3;i++)
         fprintf(stderr,"\t%f",heavyCentroid[i]);

      fprintf(stderr,"\nHeavy chain projection:\t");

      for(i=0;i<3;i++)
         fprintf(stderr,"\t%f",heavyProjection[i]);

      fprintf(stderr,"\n\n");
   }

   /* Now that the projected points have been found, find the torsion
      angle using the function "blPhi"
   */
   torsionAngle=blPhi(lightProjection[0],
                      lightProjection[1],
                      lightProjection[2],
                      lightCentroid[0],
                      lightCentroid[1],
                      lightCentroid[2],
                      heavyCentroid[0],
                      heavyCentroid[1],
                      heavyCentroid[2],
                      heavyProjection[0],
                      heavyProjection[1],
                      heavyProjection[2]);

   return(torsionAngle);
}


/************************************************************************/
/*>void Usage()
   ------------
*//**
   Print a usage message

-  19.10.16 updated for V2.0  By: ACRM
*/
void Usage()
{
   fprintf(stderr, "\nabpackingangle V2.0 (c) 2007-2016, UCL, Abhi \
Raghavan and Andrew Martin\n");
   
   fprintf(stderr, "\nUsage: abpackingangle [-p pdbcode][-o vecfile][-v]\
[-q] [in.pdb [out.txt]]\n");
   fprintf(stderr, "           -p Specify a PDB code to be printed with \
the results\n");
   fprintf(stderr, "           -o Create a PDB file containing the \
vectors used for angle calculations\n");
   fprintf(stderr, "           -v Verbose\n");
   fprintf(stderr, "           -q Quiet - prints only the angle\n");

   fprintf(stderr, "\nabpackingangle calculates the packing angle \
between VH and VL domains\n");
   fprintf(stderr, "as described by Abhinandan and Martin 23(2010),\
689-697.\n\n");
}


/************************************************************************/
/*>BOOL ParseCmdLine(int argc, char **argv, char *inFilename, 
                         char *outFilename, char *pdbCode, char *vecFile,
                         BOOL *verbose, BOOL *quiet, BOOL *printDistance,
                         BOOL *force)
   --------------------------------------------------------------------
*//**
   Parse the command line. This is a complete rewrite for V2.0

-  19.10.16 Original   By: ACRM
-  19.06.19 Added -d and -f
*/
BOOL ParseCmdLine(int argc, char **argv, char *inFile, 
                  char *outFile, char *pdbCode, char *vecFile,
                  BOOL *verbose, BOOL *quiet, BOOL *printDistance, 
                  BOOL *force)
{
   inFile[0]  = '\0';
   outFile[0] = '\0';
   pdbCode[0] = '\0';
   vecFile[0] = '\0';
   *quiet     = FALSE;
   *verbose   = FALSE;
   
   argc--;
   argv++;
   
   while(argc)
   {
      if(argv[0][0] == '-')
      {
         switch(argv[0][1])
         {
         case 'h':
            return(FALSE);
            break;
         case 'q':
            *quiet = TRUE;
            break;
         case 'p':
            argc--;
            argv++;
            if(!argc)
               return(FALSE);
            strncpy(pdbCode, argv[0], MAXPDBCODE);
            pdbCode[MAXPDBCODE-1] = '\0';
            break;
         case 'o':
            argc--;
            argv++;
            if(!argc)
               return(FALSE);
            strncpy(vecFile, argv[0], MAXFILENAME);
            vecFile[MAXPDBCODE-1] = '\0';
            break;
         case 'v':
            *verbose = TRUE;
            break;
         case 'd':
            *printDistance = TRUE;
            break;
         case 'f':
            *force = TRUE;
            break;
         default:
            return(FALSE);
         }
      }
      else
      {
         if(argc > 2)
            return(FALSE);
         
         /* Copy the first to inFile                                    */
         if(argc)
         {
            strcpy(inFile, argv[0]);
            argc--;
         }

         /* Copy the second to outFile                                  */
         if(argc)
         {
            strcpy(outFile, argv[0]);
            argc--;
         }

         return(TRUE);
      }

      argc--;
      argv++;

   }
   
   return(TRUE);
}


/************************************************************************/
/*>BOOL Plot(FILE *fpIn, FILE *fpOut, FILE *fpVecOut,
             BOOL displayOutputFlag, BOOL verbose, BOOL quiet,
             char **lightKeyResidues, char **heavyKeyResidues,
             BOOL printDistance, BOOL force)
   --------------------------------------------------------
*//**
   The main routine for doing the calculation
*/
BOOL Plot(FILE *fpIn, FILE *fpOut, FILE *fpVecOut,
          BOOL displayOutputFlag, BOOL verbose, BOOL quiet,
          char **lightKeyResidues, char **heavyKeyResidues,
          BOOL printDistance, BOOL force)
{
   int  lightChainNumberOfAtoms = 0,
        i                       = 0,
        nres;
   PDB  *pdb=NULL;

   REAL **lightKeyResidueCoords,
        **heavyKeyResidueCoords,
        lightEigenVector[3],
        heavyEigenVector[3],
        lightChainCentroid[3],
        heavyChainCentroid[3],
        **modifiedLightCoordinates = NULL,
        **modifiedHeavyCoordinates = NULL,
        torsionAngle               = 0.0,
        dist                       = 0.0;

   /* Create coordinate arrays                                          */
   if((lightKeyResidueCoords = (REAL **)blArray2D(sizeof(REAL), 
                                                  NUMRESPOS, 3))==NULL)
   {
      fprintf(stderr,"Error (abpackingangle): No memory for light \
chain key residue coordinates\n");
      return(FALSE);
   }
   if((heavyKeyResidueCoords = (REAL **)blArray2D(sizeof(REAL), 
                                                  NUMRESPOS, 3))==NULL)
   {
      fprintf(stderr,"Error (abpackingangle): No memory for heavy \
chain key residue coordinates\n");
      return(FALSE);
   }
   if((modifiedLightCoordinates = (REAL **)blArray2D(sizeof(REAL),
                                                     (NUMRESPOS/2),
                                                     3))==NULL)
   {
      fprintf(stderr,"Error (abpackingangle): No memory for light \
chain modified coordinates\n");
      return(FALSE);
   }
   if((modifiedHeavyCoordinates = (REAL **)blArray2D(sizeof(REAL),
                                                     (NUMRESPOS/2),
                                                     3))==NULL)
   {
      fprintf(stderr,"Error (abpackingangle): No memory for heavy \
chain modified coordinates\n");
      return(FALSE);
   }
   
   /* Read PDB file                                                     */
   pdb=blReadPDB(fpIn,&lightChainNumberOfAtoms);

   /* Extract coordinates of key residue CAs                            */
   if((nres =  ReadCoordinates(lightKeyResidues, lightKeyResidueCoords,
                               pdb)) != NUMRESPOS)
   {
      return(FALSE);
   }
   if((nres = ReadCoordinates(heavyKeyResidues, heavyKeyResidueCoords,
                              pdb)) != NUMRESPOS)
   {
      return(FALSE);
   }

   /* Find mid points of atoms that are structurally adjacent.
    * Light chain....
    */
   for(i=0; i<NUMRESPOS/2; i++)
   {
      modifiedLightCoordinates[i][0] = 
         (lightKeyResidueCoords[i][0] +
          lightKeyResidueCoords[i+(NUMRESPOS/2)][0])/2;

      modifiedLightCoordinates[i][1] =
         (lightKeyResidueCoords[i][1] +
          lightKeyResidueCoords[i+(NUMRESPOS/2)][1])/2;

      modifiedLightCoordinates[i][2] =
         (lightKeyResidueCoords[i][2] +
          lightKeyResidueCoords[i+(NUMRESPOS/2)][2])/2;
   }

   /* Heavy chain...                                                    */
   for(i=0; i<NUMRESPOS/2; i++)
   {
      modifiedHeavyCoordinates[i][0] = 
         (heavyKeyResidueCoords[i][0] +
          heavyKeyResidueCoords[i+(NUMRESPOS/2)][0])/2;
      
      modifiedHeavyCoordinates[i][1] =
         (heavyKeyResidueCoords[i][1] +
          heavyKeyResidueCoords[i+(NUMRESPOS/2)][1])/2;

      modifiedHeavyCoordinates[i][2] =
         (heavyKeyResidueCoords[i][2] +
          heavyKeyResidueCoords[i+(NUMRESPOS/2)][2])/2;
   }

   /* Fit the least squares distant line for the points.                */
   ComputeBestFitLine(modifiedLightCoordinates, NUMRESPOS/2, 3,
                      lightChainCentroid, lightEigenVector);

   ComputeBestFitLine(modifiedHeavyCoordinates, NUMRESPOS/2, 3,
                      heavyChainCentroid, heavyEigenVector);

   /* Calculate the torsion angle using the following points:

      Light chain: Light chain centroid and mid point of L35 and L88.
      Heavy chain: Heavy chain centroid and mid point of H36 and H92.
   */

   dist = 0.0;
   for(i=0; i<3; i++)
   {
      dist += (lightChainCentroid[i] - heavyChainCentroid[i]) *
              (lightChainCentroid[i] - heavyChainCentroid[i]);
   }
   dist = sqrt(dist);

   if((dist <= MAXDIST) || force)
   {
      if(printDistance)
         fprintf(stdout, "Distance: %7.3f   ", dist);

      torsionAngle=calculate_torsion_angle(lightEigenVector,
                                           lightChainCentroid,
                                           modifiedLightCoordinates[0],
                                           heavyEigenVector,
                                           heavyChainCentroid,
                                           modifiedHeavyCoordinates[0], 
                                           verbose);

      torsionAngle=(torsionAngle * 180)/PI;

      /* Print the torsion angle and write the coordinates of the
         imaginary regression line into a PDB file if required.
      */
      if(quiet)
      {
         fprintf(fpOut, "%f\n",torsionAngle);
      }
      else
      {
         fprintf(fpOut, "Packing angle: %f\n",torsionAngle);
      }
   
      if(displayOutputFlag)
      {
         DrawRegressionLine(lightKeyResidueCoords, lightEigenVector,
                            NUMRESPOS, "X", fpVecOut);
         DrawRegressionLine(heavyKeyResidueCoords, heavyEigenVector,
                            NUMRESPOS, "Y", fpVecOut);
         blWritePDB(fpVecOut,pdb);
      }
   }
   else
   {
      fprintf(stderr,"Error (abpackingangle): Light and Heavy chains are \
too far apart! (%.2fA)\n", dist);
   }

   /* Free memory                                                       */
   blFreeArray2D((void *)modifiedLightCoordinates, (NUMRESPOS/2), 3);
   blFreeArray2D((void *)modifiedHeavyCoordinates, (NUMRESPOS/2), 3);
   
   blFreeArray2D((void *)lightKeyResidueCoords, NUMRESPOS, 3);
   blFreeArray2D((void *)heavyKeyResidueCoords, NUMRESPOS, 3);
   
   FREELIST(pdb, PDB);

   return(TRUE);
}


/************************************************************************/
/*>int ReadCoordinates(char **keyResidueList, REAL **coordArray, PDB *pdb)
   -----------------------------------------------------------------------
*//**
 Get the X, Y, Z coordinates of the CA atoms in light and heavy
 chain constant residue list.
*/
int ReadCoordinates(char **keyResidueList, REAL **coordArray,
                     PDB *pdb)
{
   int i  = 0;
   PDB *p = NULL;

   for(i=0; i<NUMRESPOS; i++)
   {
      if((p = blFindResidueSpec(pdb, keyResidueList[i]))==NULL)
      {
         fprintf(stderr,"Error (abpackingangle): Residue not found \
(%s)\n", keyResidueList[i]);
         return(0);
      }
      
      if((p=blFindAtomInRes(p, "CA"))==NULL)
         break;

      coordArray[i][0]=(REAL)p->x;
      coordArray[i][1]=(REAL)p->y;
      coordArray[i][2]=(REAL)p->z;
   }

   return i;
}
