# include <string.h>
# include "arrays.h"

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

} /* End of function "find_largest_smallest_indices_double". */


void free_array_2D_char(char **array,int numberOfElements)
{
   int i=0;

   for(i=0;i<numberOfElements;i++)
   {
      if(array[i] != NULL)
      {
	 free(array[i]);
      }
   }

   if(array != NULL)
   {
      free(array);
   }

} /* End of function "free_array_2D_char" */


void free_array_2D_double(double **array,int numberOfElements)
{
   int i=0;

   for(i=0;i<numberOfElements;i++)
   {
      if(array[i] != NULL)
      {
	 free(array[i]);
      }
   }

   if(array != NULL)
   {
      free(array);
   }

} /* End of function "free_array_2D_char" */


void free_array_2D_int(int **array,int numberOfElements)
{
   int i=0;

   for(i=0;i<numberOfElements;i++)
   {
      if(array[i] != NULL)
      {
	 free(array[i]);
      }
   }

   if(array != NULL)
   {
      free(array);
   }

} /* End of function "free_array_2D_int" */


void free_array_2D(void **array,int numberOfElements)
{
   int i=0;

   for(i=0;i<numberOfElements;i++)
   {
      if(array[i] != NULL)
      {
	 free(array[i]);
      }
   }

   if(array != NULL)
   {
      free(array);
   }

} /* End of function "free_array_2D_char" */


void q_sort(double *numbers, int left, int right)
{
   double pivot, l_hold, r_hold;
 
   l_hold = left;
 
   r_hold = right;
   
   pivot = numbers[left];
   
   while (left < right)
   {
      while ((numbers[right] >= pivot) && (left < right))
      {
	 right--;
      }

      if (left != right)
      {
	 numbers[left] = numbers[right];
	 left++;
      }

      while ((numbers[left] <= pivot) && (left < right))
      {
	 left++;
      }

      if (left != right)
      {
	 numbers[right] = numbers[left];
	 right--;
      }
   }
   
   numbers[left] = pivot;
   pivot = left;
   left = l_hold;
   right = r_hold;
   
   if (left < pivot)
   {
      q_sort(numbers, left, pivot-1);
   }
   if (right > pivot)
   {
      q_sort(numbers, pivot+1, right);
   }

} /* End of function "q_sort" */


/* Implementation of quick sort to sort an array. Code reproduced from:

   http://linux.wku.edu/~lamonml/algor/sort/quick.html
*/

void quick_sort(double *numbers, int array_size)
{
   q_sort(numbers, 0, array_size - 1);

} /* End of function "quick_sort" */



void initialise_double_array(double *array, int numberOfElements, double value)
{
   int i=0;

   for(i=0;i < numberOfElements; i++)
   {
      array[i] = value;
   }

} /* End of function "initialise_double_array" */


void initialise_int_array(int *array, int numberOfElements, int value)
{
   int i=0;

   for(i=0;i < numberOfElements; i++)
   {
      array[i] = value;
   }

} /* End of function "initialise_int_array" */


void initialise_array_of_strings(char **array, int numberOfStrings, char *stringToInitialise)
{
   int i = 0;

   for(i=0; i < numberOfStrings; i++)
   {
      strcpy(array[i],stringToInitialise);
   }

} /* End of function "initialise_array_of_strings" */


void convert_to_upper_case(char *sourceString, char *destString)
{
   int i = strlen(sourceString);

   for(i = 0; i < strlen(sourceString); i++)
   {
      destString[i] = toupper(sourceString[i]);
   }

   destString[i] = '\0';

} /* End of function "convert_to_upper_case" */


void convert_to_lower_case(char *sourceString, char *destString)
{
   int i = 0;

   for(i = 0; i < strlen(sourceString); i++)
   {
      destString[i] = tolower(sourceString[i]);
   }

   destString[i] = '\0';

} /* End of function "convert_to_lower_case" */
