# include <stdio.h>
# include <stdlib.h>
# include "bioplib/MathType.h"
# include "standard.h"

#ifndef ARRAYS
#define ARRAYS

#ifndef FREE
# define FREE(X)\
   do\
   {\
      if( (X) != NULL )\
      {\
	 free((X));\
	 (X)=NULL;\
      }\
\
   }while(0);
#endif

void find_largest_smallest_indices_double(double *doubleArray,
                                          int numberOfElements,
                                          int *smallestValueIndex,
                                          int *largestValueIndex);

void free_array_2D_char(char **array,int numberOfElements);

void free_array_2D_double(double **array,int numberOfElements);

void free_array_2D_int(int **array,int numberOfElements);

void free_array_2D(void **array,int numberOfElements);

void quick_sort(double *numbers, int array_size);

void q_sort(double *numbers, int left, int right);

void initialise_double_array(double *array, int numberOfElements, double value);

void initialise_int_array(int *array, int numberOfElements, int value);

void initialise_array_of_strings(char **array, int numberOfStrings, char *stringToInitialise);

void convert_to_upper_case(char *sourceString, char *destinationString);

void convert_to_lower_case(char *sourceString, char *destinationString);

#endif
