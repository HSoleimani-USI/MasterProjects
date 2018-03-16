/* ==================================================================== *
 *								        *
 *  block_dgemm.c -- Implemant a block matrix multiplication routine    *
 *                                                                      *
 * ==================================================================== */

//  Author: Hanieh Soleimani 


#include "square_dgemm.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* block parameter ... */
#ifndef BLOCK_SIZE
#  define BLOCK_SIZE ((unsigned) 16)
#endif


#define min(a,b) (((a) < (b)) ? (a) : (b))




/**
 *  square_dgemm -- multiply two block matrices A and B adding result to C, result is C = C + A*B
 */


void square_dgemm (const double  *A, const double  *B,  double  *C, const unsigned  M)
{

	/* TODO: implement the blocked matrix-matrix multiplication */
	   // initializing the result matrix c
	   // for(int i=0; i<=M; ++i)
	   // 	  for(int j=0; i<=M; ++j)
	   // 	  	 C[i][j]=0;
	  int i, j,k;
	  int inner_i,inner_j,inner_k;

	    for(int k=0; k<=M; k+=BLOCK_SIZE)
	    	for(int j=0; j<=M; j+=BLOCK_SIZE)
	    		for(int i=0; i<=M; i+=BLOCK_SIZE)

	    			for(int inner_k=k; inner_k< min(k + BLOCK_SIZE, M); ++inner_k)
	    				for (int inner_j =j; inner_j < min(j + BLOCK_SIZE, M); ++inner_j)
	    					for (int inner_i = i; inner_i < min(i + BLOCK_SIZE, M); ++inner_i)
	    					   C[inner_k*M+inner_j] += A[inner_k*M+inner_i] * B[inner_i*M+inner_j];


        // printf("implement square_dgemm routine\n");
        // exit(0);
}

