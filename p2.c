/*
	Name: Kirolous shihataa
	instructor: Ralph butler  
	Due date: 
	Assignment: p2 
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#include <pthread.h>
#include <sys/time.h>

int numRows, numCols ; 
float *A[5000], *B[5000];
float epsilon;
float maxDiffs[64]; 
int nelems_per_thread;
pthread_barrier_t barr; 
int nthreads; 
void* sub(void *arg);

int main( int argc, char *argv[])
{
	if(argc != 9)
	{
		printf("incorrect number of arguments, given %d arguments \n", argc);
		return 0; 
	}

	numRows = atoi(argv[1]); 
	numCols = atoi(argv[2]); 
	float topTemp = atof(argv[3]); 
	float leftTemp = atof(argv[4]); 
	float rightTemp = atof(argv[5]); 
	float botTemp = atof(argv[6]); 
	epsilon  = atof(argv[7]); 
	nthreads  = atoi(argv[8]);
	int ids[64]; 
	float *T[numRows];

	struct timeval t0, t1, dt; 
	pthread_t threadid[64];

	int i, j, x, y; 
	for(i = 0; i < numRows; i++ )
	{
		A[i] = (float *)malloc((numCols) * sizeof(float)); 
		B[i] = (float *)malloc((numCols) * sizeof(float)); 
		T[i] = (float *)malloc((numCols) * sizeof(float)); 
	}
	for(i=0; i < numRows; i++)
	{
		for (j=0; j<numCols; j++)
		{
			A[i][j] = 0.0;
			B[i][j] = 0.0; 
		}

	}
	for (j = 0 ; j < numCols; j++  )
	{	
		A[0][j] = topTemp; 
		B[0][j] = topTemp; 
	}
	for (i  = 0 ; i < numRows; i++ )
	{	
		A[i][0] = leftTemp; 
		B[i][0] = leftTemp; 
	}
	for (i = 0 ; i < numRows; i++ )
	{
		A[i][numCols-1] = rightTemp; 
		B[i][numCols-1] = rightTemp; 
	}
	for (j = 0 ; j < numCols; j++ )
	{
		A[numRows-1][j] = botTemp; 
		B[numRows-1][j] = botTemp; 
	}
        
  
	float edgeSum = 0.0;
	int numEdgeVals = 0;
	
	for(j = 0; j < numCols; j++)
	{
		edgeSum += A[0][j] + A[numRows-1][j];
		numEdgeVals += 2;
	}
	for (i = 1; i < numRows-1; i++ )
	{
		edgeSum+= A[i][0] + A[i][numCols-1]; 
		numEdgeVals += 2;
	}

	float interiorInitVal = edgeSum / numEdgeVals; 
	printf("%15.6f <---- edge sumd (DBG) \n", edgeSum);
	printf("%15.6f <---- inital val for all interior points (DBG) \n", interiorInitVal );

	int rowidx, epoch, colidx;

	for (rowidx = 1 ; rowidx < numRows-1 ; rowidx ++)
		for (colidx = 1 ; colidx < numCols-1 ;colidx ++)
			A[rowidx][colidx] = interiorInitVal; 
	

	// number of elements per thread
	
	nelems_per_thread = (numRows / nthreads); 



	pthread_barrier_init(&barr, NULL, nthreads ); // init barrier		
	
	gettimeofday(&t0, 0);
	for (i = 1 ; i < nthreads; i++){
		ids[i] = i;
		pthread_create(&threadid[i], NULL, sub, &ids[i]);
	}

	// main thread work 
	float maxDiff; 
	float sumNeighbors; 
	float newVal;
       	float absDiff;
	
	int begin = (0 * nelems_per_thread) + 1;
	int end ; 
	
	if(0 == nthreads-1 )
		end  = numRows - 1; 
	else
		end = begin + nelems_per_thread; 
		

//	printf ("ID0 beginning %d  ending %d \n",begin,  end);	
	for(epoch =0; epoch<5000; epoch++)
	{
		maxDiff = 0.0; 
		for (rowidx = begin; rowidx < end; rowidx++)
		{
			for (colidx = 1; colidx < (numCols-1); colidx++)
			{
				sumNeighbors = A[rowidx-1][colidx] + A[rowidx+1][colidx] + A [rowidx][colidx-1] + A[rowidx][colidx+1]; 
				newVal = sumNeighbors / 4.0; 
				absDiff = fabsf(newVal - A[rowidx][colidx]) ; 
				
				if(absDiff > maxDiff)
					maxDiff = absDiff; 
				B[rowidx][colidx] = newVal;
			}
		}

	
		// wait for the children to do the work 
		pthread_barrier_wait(&barr);
	
		float logE = log(epoch)/log(2); 
		if(epoch && ceilf(logE) == logE )
			printf("%3d %8.6f \n", epoch, maxDiff);

		for(x = 1; x < numRows-1; x++)
		{
				T[x]= B[x];
				B[x]= A[x];
				A[x]= T[x]; 
		}
		

		// find which thread had the max diff
		for(i = 1 ; i < nthreads; i++)
		{
			if(maxDiffs[i]> maxDiff)
				maxDiff = maxDiffs[i];
		}
		
		maxDiffs[0] = maxDiff;
		pthread_barrier_wait(&barr); // let the threads hang till 

		if(maxDiff < epsilon)
		{
			printf("%3d %8.6f \n", epoch, maxDiff);
			break; 
		}
		

	}

	for(i = 1 ; i < nthreads ; i++)
		pthread_join(threadid[i],NULL);
	gettimeofday(&t1, 0);
printf("\nTOTAL TIME  %.2f seconds\n\n",
      (((t1.tv_sec * 1000000 + t1.tv_usec) -
    (t0.tv_sec * 1000000 + t0.tv_usec))*1.0)/1000000.0);	
return 0;
}

void *sub(void *arg)
{
	float maxDiff; 
	float sumNeighbors; 
	float newVal;
       	float absDiff;
	
	int myid = *((int *)arg);
	int begin = (myid * nelems_per_thread) + 1;
	int end ; 
	int epoch; 
	int rowidx, colidx;

	if(myid == nthreads-1 )
		end  = numRows - 1; 
	else
		end = begin + nelems_per_thread; 
		
//	printf ("ID%d beginning %d  ending %d \n",myid, begin,  end);
	for(epoch = 0; epoch<5000; epoch++)
	{
		maxDiff = 0.0; 
		for (rowidx = begin; rowidx < end ; rowidx++)
		{
			for (colidx = 1; colidx < (numCols-1); colidx++)
			{
				sumNeighbors = A[rowidx-1][colidx] + A[rowidx+1][colidx] + A [rowidx][colidx-1] + A[rowidx][colidx+1]; 
				newVal = sumNeighbors / 4.0; 
				absDiff = fabsf(newVal - A[rowidx][colidx]) ; 
				
				if(absDiff > maxDiff)
					maxDiff = absDiff; 
				B[rowidx][colidx] = newVal;
			}
		}

		maxDiffs[myid] = maxDiff;
		pthread_barrier_wait(&barr);  // join the barrier and wait for everybody to do the work
		pthread_barrier_wait(&barr); // join the barrier and wait for the main thread to do the sawp and to calc maxdiff 

		if(maxDiffs[0] < epsilon)
			break; 	
	} 
}
