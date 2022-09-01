//gcc -O3 -lm -lrt -lfftw3 -o STRUCTURE_FACTOR.EXE test.c &&  ./STRUCTURE_FACTOR.EXE

#include <stdio.h>  
#include <stdlib.h>
#include <string.h>
#include <math.h>


#define TPI 6.283185307179586477

#define WIDTH	5
/* Create abstract data type for vector */
typedef struct {
	double q[4];
	double p_s;
	double p_m;
} pot_q;

// A utility function to swap two elements 
void swap(pot_q* a, pot_q* b) 
{ 
    pot_q t = *a; 
    *a = *b; 
    *b = t; 
} 
  
/* This function takes last element as pivot, places 
   the pivot element at its correct position in sorted 
    array, and places all smaller (smaller than pivot) 
   to left of pivot and all greater elements to right 
   of pivot */
int partition (pot_q arr[], int low, int high) 
{ 
    pot_q pivot = arr[high];    // pivot 
    int i = (low - 1);  // Index of smaller element 
  
    for (int j = low; j <= high- 1; j++) 
    { 
        // If current element is smaller than the pivot 
        if (arr[j].q[3] < pivot.q[3]) 
        { 
            i++;    // increment index of smaller element 
            swap(&arr[i], &arr[j]); 
        } 
    } 
    swap(&arr[i + 1], &arr[high]); 
    return (i + 1); 
} 
  
/* The main function that implements QuickSort 
 arr[] --> Array to be sorted, 
  low  --> Starting index, 
  high  --> Ending index */
void quickSort(pot_q arr[], int low, int high) 
{ 
    if (low < high) 
    { 
        /* pi is partitioning index, arr[p] is now 
           at right place */
        int pi = partition(arr, low, high); 
  
        // Separately sort elements before 
        // partition and after partition 
        quickSort(arr, low, pi - 1); 
        quickSort(arr, pi + 1, high); 
    } 
} 








int main()  
{  
    int 			i;  
                     
    pot_q			a[100000]; 
      
    	
	for(i = 0; i < 100000; i++){
		a[i].q[3] = rand()/100000.0; 
		a[i].p_s  = i;
		a[i].p_m  = i;
		}
		
	quickSort(a,0,99999);
	for(i = 0; i < 100; i++){
		printf ("%12.8lf,	", a[i].q[3]) ;
		printf ("%12.8lf,	", a[i].p_s );
		printf ("%12.8lf\n", a[i].p_m);
	}
	printf("\n**************************************\n\n");
	


    return 0; 
 }
 
 


  
