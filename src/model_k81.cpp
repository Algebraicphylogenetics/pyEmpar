/*
 *  model_k81.cpp
 *  
 *  Created by Ania M. Kedzierska on 11/11/11.
 *  Copyright 2011 Politecnic University of Catalonia, Center for Genomic Regulation.  This is program can be redistributed, modified or else as given by the terms of the GNU General Public License. 
 *  
 */

#include "model_k81.h"
#include "random.h"
#include <cstdlib>
#include <iostream>
#include <cmath>
#include <algorithm>

#include "miscelania.h"


//////////////////////////////////////////////////////////////////
// This file implements the functions specific to K81*, 4 states
//////////////////////////////////////////////////////////////////


// Fills in a K81 matrix.
void K81_matrix(double a, double b, double c, double d, TMatrix &tm) {
  tm[0][0] = a;
  tm[0][1] = b;
  tm[0][2] = c;
  tm[0][3] = d;

  tm[1][0] = b;
  tm[1][1] = a;
  tm[1][2] = d;
  tm[1][3] = c;

  tm[2][0] = c;
  tm[2][1] = d;
  tm[2][2] = a;
  tm[2][3] = b;

  tm[3][0] = d;
  tm[3][1] = c;
  tm[3][2] = b;
  tm[3][3] = a;
}


long K81_matrix_structure(long i, long j) {
  if (i == j) return 0;
  else if ((i+j) % 4 == 1) return 1;
  else if ((i+j) % 2 == 0) return 2;
  else if (i+j == 3) return 3;
  else {
    std::cout << "ERROR: unknown condition in K81_matrix_structure." << std::endl;
    exit(1);
  }
}


long K81_root_structure(long i) {
  return 0;
}


Permutation K81_perm(long a, long b, long c, long d) {
  Permutation p;
  p.resize(4);
  p[0] = a; p[1] = b; p[2] = c; p[3] = d;
  return p;
}

void K81_list_permutations(std::list<Permutation> &L) {
  L.clear();
  L.push_back(K81_perm(0,1,2,3));   // det 1
  L.push_back(K81_perm(1,0,3,2));   // det 1
  L.push_back(K81_perm(2,3,0,1));   // det 1
  L.push_back(K81_perm(3,2,1,0));   // det 1
}




// Computes the mle for the transition matrix of a given edge from the
// the marginalization of the counts on that edge.
// N is a matrix representing the counts on that edge
// tm is the matrix where the mle is stored. tm must be a stochastic matrix.
void K81_mle_edge(TMatrix &N, TMatrix &tm) {
  double a, b, c, d, sum;

  a = N[0][0] + N[1][1] + N[2][2] + N[3][3];
  b = N[0][1] + N[1][0] + N[2][3] + N[3][2];
  c = N[0][2] + N[2][0] + N[1][3] + N[3][1];
  d = N[0][3] + N[3][0] + N[1][2] + N[2][1];
  sum = a+b+c+d;
  a = a/sum;
  b = b/sum;
  c = c/sum;
  d = d/sum;

  K81_matrix(a,b,c,d,tm);
}

// Computes the mle for the root distribution of a given edge from the
// marginalization of the counts to the root node.
// N is a vector representing the counts on the root
// r is a vector where the root distribution is stored.
void K81_mle_root(Root &s, Root &r) {
  long i;
  // For K81 the root distribution is uniform.
  for (i=0; i < 4; i++) {
    r[i] = 0.25;
  }
}


//K81 root distribution:
void K81_random_root(Root &r) {
  r[0] = 0.25;
  r[1] = 0.25;
  r[2] = 0.25;
  r[3] = 0.25;
}

// Produces a random K81 transition matrix with diagonal entries >= 0.7.
void K81_random_edge(TMatrix &tm) {
  double a,b,c,d;
  double dmin = 0.7;

  do {
    a = uniform_real(dmin, 1);
    b = uniform_real(0, 1-dmin);
    c = uniform_real(0, 1-dmin);
    d = 1 - a - b - c;
  } while (d < 0);


  K81_matrix(a,b,c,d,tm);
}


// Random K81 transition matrix for a given branch length
void K81_random_edge_length(double len, TMatrix &tm) {
  double a,b,c,d;
  double K, x1, x2, x3, s, m, M, p, q;
  K = exp(-4*len);

  // s is the only real root of x(x+1)^2 - 4K = 0
  s = - 2./3. + 1./3.*pow(1. + 54.*K + 6.*sqrt(3.*K + 81.*K*K), 1./3.)
              + 1./3.*pow(1. + 54.*K - 6.*sqrt(3.*K + 81.*K*K), 1./3.);

  if (s > 1 || s < 0) {
    std::cout << "Error: Something strange happened at K81_random_edge_length. Bad s.";
    std::cout << std::endl;
  }
  x1 = uniform_real(s, 1);
  p = sqrt((x1*(x1-1)*(x1-1)+4*K)/x1);
  q = sqrt((x1*(x1+1)*(x1+1)-4*K)/x1);
  m = 0.5*(x1 + std::max<double>(1-q, p-1));
  M = 0.5*(1 + std::min<double>(p-x1, q+x1));

  if (m > M || m < 0) {
    std::cout << "Error: Something strange happened at K81_random_edge_length: bad bounds m and M.";
    std::cout << std::endl;
  }
  x2 = uniform_real(m, M);

  double gamma;
  gamma = uniform_real(0,1);
  if (gamma > 0.5) x1 = -x1;

  gamma = uniform_real(0,1);
  if (gamma > 0.5) x2 = -x2;

  x3 = K/(x1*x2);
  b = 0.25*(1 - x1 - x2 + x3);
  c = 0.25*(1 - x1 + x2 - x3);
  d = 0.25*(1 + x1 - x2 - x3);
  a = 1 - b - c - d;

  K81_matrix(a,b,c,d,tm);
}


// Random biologically meaningful K81 transition matrix with given length.
// Here biologically meaningful means that diagonal entries are maximal in the column they belong to (and the row in this case)
// ( Chang's DLC condition).
void K81_random_edge_bio_length(double len, TMatrix &tm) {
  TMatrix tmaux;
  long i0, i;

  tmaux.resize(4);
  for(i=0; i < 4; i++) {
    tmaux[i].resize(4);
  }

  K81_random_edge_length(len, tmaux);

  // Permute with the row that starts with the largest entry in column.
  // All 4 possible permutations have det = 1
  i0 = max_in_col(tmaux, 0);
  K81_matrix(tmaux[i0][0], tmaux[i0][1], tmaux[i0][2], tmaux[i0][3], tm);
}
