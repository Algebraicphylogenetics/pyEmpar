/*
 *  seqUtil.h
 *  
 *  Created by Ania M. Kedzierska on 11/11/11.
 *  Copyright 2011 Politecnic University of Catalonia, Center for Genomic Regulation.  This is program can be redistributed, modified or else as given by the terms of the GNU General Public License. 
 *  
 */

#ifndef __SEQUTIL_H__
#define __SEQUTIL_H__

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define APPEND_LEN	256

#ifdef __SEQUTIL_C__
	void seqMemInit();
	void* seqMalloc(int size);
	void seqFree();
	void seqFreeAll();
	void inputString(char *input, char **ppcStr, int *iLen, int *iMaxLen);
#else
	extern void seqMemInit();
	extern void* seqMalloc(int size);
	extern void seqFreeAll();
	extern void seqFree();
	extern void inputString(char *input, char **ppcStr, int *iLen, int *iMaxLen);
#endif



#endif

