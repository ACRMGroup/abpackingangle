# include <stdio.h>
# include <stdlib.h>
# include "bioplib/MathType.h"

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

void find_largest_smallest_indices_REAL(REAL *REALArray,
                                          int numberOfElements,
                                          int *smallestValueIndex,
                                          int *largestValueIndex);

void free_array_2D_char(char **array,int numberOfElements);

void free_array_2D_REAL(REAL **array,int numberOfElements);

void free_array_2D_int(int **array,int numberOfElements);

void free_array_2D(void **array,int numberOfElements);

void quick_sort(REAL *numbers, int array_size);

void q_sort(REAL *numbers, int left, int right);

void initialise_REAL_array(REAL *array, int numberOfElements, REAL value);

void initialise_int_array(int *array, int numberOfElements, int value);

void initialise_array_of_strings(char **array, int numberOfStrings, char *stringToInitialise);

void convert_to_upper_case(char *sourceString, char *destinationString);

void convert_to_lower_case(char *sourceString, char *destinationString);

#endif
