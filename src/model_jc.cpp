/*
 *  model_jc.cpp
 *  
 *  Created by Ania M. Kedzierska on 11/11/11.
 *  Copyright 2011 Politecnic University of Catalonia, Center for Genomic Regulation.  This is program can be redistributed, modified or else as given by the terms of the GNU General Public License. 
 *  
 */


#include "model_jc.h"
#include "random.h"
#include <cstdlib>
#include <iostream>
#include <cmath>
#include <algorithm>



//////////////////////////////////////////////////////////////////
// This file implements specific functions for JC*, 4 states
//////////////////////////////////////////////////////////////////

// Fills in a JC matrix.
void JC_matrix(double a, double b, TMatrix &tm) {
  tm[0][0] = a;
  tm[0][1] = b;
  tm[0][2] = b;
  tm[0][3] = b;

  tm[1][0] = b;
  tm[1][1] = a;
  tm[1][2] = b;
  tm[1][3] = b;

  tm[2][0] = b;
  tm[2][1] = b;
  tm[2][2] = a;
  tm[2][3] = b;

  tm[3][0] = b;
  tm[3][1] = b;
  tm[3][2] = b;
  tm[3][3] = a;
}


long JC_matrix_structure(long i, long j) {
  if (i==j) return 0;
  else return 1;
}

long JC_root_structure(long i) {
  return 0;
}

Permutation JC_perm(long a, long b, long c, long d) {
  Permutation p;
  p.resize(4);
  p[0] = a; p[1] = b; p[2] = c; p[3] = d;
  return p;
}

void JC_list_permutations(std::list<Permutation> &L) {
  L.clear();
  L.push_back(JC_perm(0,1,2,3));
}




// Computes the mle for the transition matrix of a given edge from the
// the marginalization of the counts on that edge.
// N is a matrix representing the counts on that edge
// tm is the matrix where the mle is stored. tm must be a stochastic matrix
// of the given model !
void JC_mle_edge(TMatrix &N,  TMatrix &tm) {
  double a, b, sum;

  a = N[0][0] + N[1][1] + N[2][2] + N[3][3];
  b = N[0][1] + N[1][0] + N[2][3] + N[3][2] +
      N[0][2] + N[2][0] + N[1][3] + N[3][1] +
      N[0][3] + N[3][0] + N[1][2] + N[2][1];
  sum = a+b;
  a = a/sum;
  b = b/(3.*sum);

  JC_matrix(a,b,tm);
}

// Computes the mle for the root distribution of a given edge from the
// marginalization of the counts to the root node.
// N is a vector representing the counts on the root
// r is a vector where the root distribution is stored.
void JC_mle_root(Root &s, Root &r) {
  long i;
  // For K81 the root distribution is uniform.
  for (i=0; i < 4; i++) {
    r[i] = 0.25;
  }
}


// Produces a random JC root distribution. This is a bit silly, as the only
// possible JC root distribution is the uniform one ...
void JC_random_root(Root &r) {
  r[0] = 0.25;
  r[1] = 0.25;
  r[2] = 0.25;
  r[3] = 0.25;
}

// Produces a random JC transition matrix with diagonal entries >= 0.7.
void JC_random_edge(TMatrix &tm) {
  double a,b;
  double dmin = 0.7;

  a = uniform_real(dmin,1);
  b = (1-a)/3;

  JC_matrix(a,b,tm);
}


// Random JC transition matrix of given length.
void JC_random_edge_length(double len, TMatrix &tm) {
  double K;
  double a,b;
  K = exp(-4*len);

  b = 0.25*(1-pow(K, 1./3.));
  a = (1-3*b);

  JC_matrix(a,b,tm);
}


// For JC there's only one matrix with given length.
void JC_random_edge_bio_length(double len, TMatrix &tm) {
  JC_random_edge_length(len, tm);
}
