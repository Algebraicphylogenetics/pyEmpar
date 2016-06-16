/*
 *  model_k80.cpp
 *  
 *  Created by Ania M. Kedzierska on 11/11/11.
 *  Copyright 2011 Politecnic University of Catalonia, Center for Genomic Regulation.  This is program can be redistributed, modified or else as given by the terms of the GNU General Public License. 
 *  
 */

#include "model_k80.h"
#include "random.h"
#include <cstdlib>
#include <iostream>
#include <cmath>
#include <algorithm>

#include "miscelania.h"

//////////////////////////////////////////////////////////////////
// This file implements functions specific to K80*, 4 states
//////////////////////////////////////////////////////////////////


// Fills in a K80 matrix.
void K80_matrix(double a, double b, double c, TMatrix &tm) {
  tm[0][0] = a;
  tm[0][1] = b;
  tm[0][2] = c;
  tm[0][3] = b;

  tm[1][0] = b;
  tm[1][1] = a;
  tm[1][2] = b;
  tm[1][3] = c;

  tm[2][0] = c;
  tm[2][1] = b;
  tm[2][2] = a;
  tm[2][3] = b;

  tm[3][0] = b;
  tm[3][1] = c;
  tm[3][2] = b;
  tm[3][3] = a;
}

long K80_matrix_structure(long i, long j) {
  if (i == j) return 0;
  else if ((i+j) % 4 == 1) return 1;
  else if ((i+j) % 2 == 0) return 2;
  else if (i+j == 3) return 1;
  else {
    std::cout << "ERROR: unknown condition in K80_matrix_structure." << std::endl;
    exit(1);
  }
}

long K80_root_structure(long i) {
  return 0;
}


Permutation K80_perm(long a, long b, long c, long d) {
  Permutation p;
  p.resize(4);
  p[0] = a; p[1] = b; p[2] = c; p[3] = d;
  return p;
}

void K80_list_permutations(std::list<Permutation> &L) {
  L.clear();
  L.push_back(K80_perm(0,1,2,3));
  L.push_back(K80_perm(2,3,0,1));
}




// Computes the mle for the entries of te transition matrix of a given edge from the
// the marginalization of the counts on that edge.
// N is a matrix representing the counts on that edge
// tm is the matrix where the mle is stored. tm must be a stochastic matrix
void K80_mle_edge(TMatrix &N, TMatrix &tm) {
  double a, b, c, sum;

  a = N[0][0] + N[1][1] + N[2][2] + N[3][3];
  b = N[0][1] + N[1][0] + N[2][3] + N[3][2] +
      N[0][3] + N[3][0] + N[1][2] + N[2][1];
  c = N[0][2] + N[2][0] + N[1][3] + N[3][1];
  sum = a+b+c;
  a = a/sum;
  b = b/(2.*sum);
  c = c/sum;

  K80_matrix(a,b,c,tm);
}

// Computes the mle for the root distribution of a given edge from the
// marginalization of the counts to the root node.
// N is a vector representing the counts on the root
// r is a vector where the root distribution is stored.

// the root distribution is uniform:

void K80_mle_root(Root &s, Root &r) {
  long i;
  for (i=0; i < 4; i++) {
    r[i] = 0.25;
  }
}

void K80_random_root(Root &r) {
  r[0] = 0.25;
  r[1] = 0.25;
  r[2] = 0.25;
  r[3] = 0.25;
}

// Produces a random K80 transition matrix with diagonal entries >= 0.7.
void K80_random_edge(TMatrix &tm) {
  double a,b,c;
  double dmin = 0.7;

  do {
    a = uniform_real(dmin, 1);
    b = uniform_real(0, (1-dmin)/2);
    c = 1 - 2*b - a;
  } while (c < 0);

  K80_matrix(a,b,c,tm);
}

// Random K80 transition matrix of a given branch length.
void K80_random_edge_length(double len, TMatrix &tm) {
  double a, b, c;
  double K, x1, x2, t, s;
  K = exp(-4*len);

  // s is the only real root -2x^3+x^2+K = 0.
  s = 1./6. + 1./6.*pow(1. + 54.*K + 6.*sqrt(3.*K + 81.*K*K), 1./3.)
            + 1./6.*pow(1. + 54.*K - 6.*sqrt(3.*K + 81.*K*K), 1./3.);

  t = sqrt(K);

  if (t > s) {
    std::cout << "Error: Something strange happened at K80_random_edge_length. Bad s.";
    std::cout << std::endl;
  }
  x1 = uniform_real(t, s);

  double gamma = uniform_real(0,1);
  if (gamma > 0.5) x1 = -x1;

  x2 = K / (x1*x1);

  b = (1 - x2)/4.;
  c = (1 + x2 - 2*x1)/4.;
  a = 1 - 2*b - c;

  K80_matrix(a,b,c,tm);
}



// Random biologically meaningful K80 transition matrix of given length (acc. to Chang's DLC).
void K80_random_edge_bio_length(double len, TMatrix &tm) {
  TMatrix tmaux;
  long i0, i;

  tmaux.resize(4);
  for(i=0; i < 4; i++) {
    tmaux[i].resize(4);
  }

  // Loop until the permutation that puts the matrix into DLC preserves K80 structure.
  long timeout = 0;
  do {
    K80_random_edge_length(len, tmaux);
    i0 = max_in_col(tmaux, 0);
    timeout++;
  } while(!(i0 == 0 || i0 == 2) && timeout < 1000);

  if (timeout >= 1000) {
    std::cout << "ERROR: In sampling for K80 model. Can't generate DLC matrix of length " << len;
    std::cout << std::endl;
    exit(1);
  }

  K80_matrix(tmaux[i0][0], tmaux[i0][1], tmaux[i0][2], tm);
}
