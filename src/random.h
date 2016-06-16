/*
 *  random.h
 *  
 *  Created by Ania M. Kedzierska on 11/11/11.
 *  Copyright 2011 Politecnic University of Catalonia, Center for Genomic Regulation.  This is program can be redistributed, modified or else as given by the terms of the GNU General Public License. 
 *  
 */


#ifndef __RANDOM_H__
#define __RANDOM_H__

#include<vector>

void random_initialize();                     // Initialize with a seed derived from time
void random_initialize(unsigned int seed);    // Initialize with a given seed
double uniform_real(double a, double b);      // Uniform distribution between a and b
long uniform_int(long a, long b);             // Uniform integer discrete distribution between a and b
long discrete(std::vector<double> &p);        // Discrete distribution with given prob

#endif
