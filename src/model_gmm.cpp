/*
 *  model_gmm.cpp
 *  
 *  Created by Ania M. Kedzierska on 11/11/11.
 *  Copyright 2011 Politecnic University of Catalonia, Center for Genomic Regulation.  This is program can be redistributed, modified or else as given by the terms of the GNU General Public License. 
 *  
 */


#include "model_gmm.h"
#include "model_ssm.h"
#include "random.h"
#include "parameters.h"
#include <cstdlib>
#include <iostream>
#include <cmath>
#include <algorithm>

#include <boost/math/special_functions/erf.hpp>
#include <boost/numeric/ublas/matrix.hpp>

#include "expm.hpp"

#include "miscelania.h"

//////////////////////////////////////////////////////////////////
// This file implements functions specific to the GMM model
//////////////////////////////////////////////////////////////////

// Fills in a GMM row
void GMM_row(double a, double b, double c, double d, long r, TMatrix &tm) {
  tm[r][0] = a;
  tm[r][1] = b;
  tm[r][2] = c;
  tm[r][3] = d;
}



long GMM_matrix_structure(long i, long j) {
  return 4*i + j;
}


long GMM_root_structure(long i) {
  return i;
}

Permutation GMM_perm(long a, long b, long c, long d) {
  Permutation p;
  p.resize(4);
  p[0] = a; p[1] = b; p[2] = c; p[3] = d;
  return p;
}

void GMM_list_permutations(std::list<Permutation> &L) {
  L.clear();
  L.push_back(GMM_perm(0,1,2,3));
  L.push_back(GMM_perm(0,1,3,2));
  L.push_back(GMM_perm(0,2,1,3));
  L.push_back(GMM_perm(0,2,3,1));
  L.push_back(GMM_perm(0,3,1,2));
  L.push_back(GMM_perm(0,3,2,1));

  L.push_back(GMM_perm(1,0,2,3));
  L.push_back(GMM_perm(1,0,3,2));
  L.push_back(GMM_perm(1,2,0,3));
  L.push_back(GMM_perm(1,2,3,0));
  L.push_back(GMM_perm(1,3,0,2));
  L.push_back(GMM_perm(1,3,2,0));

  L.push_back(GMM_perm(2,0,1,3));
  L.push_back(GMM_perm(2,0,3,1));
  L.push_back(GMM_perm(2,1,0,3));
  L.push_back(GMM_perm(2,1,3,0));
  L.push_back(GMM_perm(2,3,0,1));
  L.push_back(GMM_perm(2,3,1,0));

  L.push_back(GMM_perm(3,0,1,2));
  L.push_back(GMM_perm(3,0,2,1));
  L.push_back(GMM_perm(3,1,0,2));
  L.push_back(GMM_perm(3,1,2,0));
  L.push_back(GMM_perm(3,2,0,1));
  L.push_back(GMM_perm(3,2,1,0));
}



// Computes the mle for the transition matrix of a given edge from the
// the marginalization of the counts on that edge.
// N is a matrix representing the counts on that edge
// tm is the matrix where the mle is stored. tm must be a stochastic matrix
void GMM_mle_edge(TMatrix &N, TMatrix &tm) {
    double s[4];
    long i,j;

    for (i=0; i < 4; i++) {
      s[i] = 0;
      for(j=0; j < 4; j++) {
        s[i] = s[i] + N[i][j];
      }
    }

    tm[0][0] = N[0][0] / s[0];
    tm[0][1] = N[0][1] / s[0];
    tm[0][2] = N[0][2] / s[0];
    tm[0][3] = N[0][3] / s[0];

    tm[1][0] = N[1][0] / s[1];
    tm[1][1] = N[1][1] / s[1];
    tm[1][2] = N[1][2] / s[1];
    tm[1][3] = N[1][3] / s[1];

    tm[2][0] = N[2][0] / s[2];
    tm[2][1] = N[2][1] / s[2];
    tm[2][2] = N[2][2] / s[2];
    tm[2][3] = N[2][3] / s[2];

    tm[3][0] = N[3][0] / s[3];
    tm[3][1] = N[3][1] / s[3];
    tm[3][2] = N[3][2] / s[3];
    tm[3][3] = N[3][3] / s[3];
}

// Computes the mle for the root distribution of a given edge from the
// marginalization of the counts to the root node.
// s is a vector representing the counts on the root
// r is a vector where the root distribution is stored.
void GMM_mle_root(Root &s, Root &r) {
  long i;
  double sum;
  sum = s[0] + s[1] + s[2] + s[3];

  for (i=0; i < 4; i++) {
    r[i] = s[i] / sum;
  }
}

// Uniform random stochastic vector of length 3.
void GMM_random_stochastic_vector3(std::vector<double> &v) {
  v[1] = uniform_real(0,1);
  v[2] = uniform_real(0,1);

  if (v[1] + v[2] > 1) {
    v[1] = 1 - v[1];
    v[2] = 1 - v[2];
  }

  v[0] = 1 - v[1] - v[2];
}

// Uniform random stochastic vector of length 4.
void GMM_random_stochastic_vector4(std::vector<double> &v) {
  double s;

  // A uniform point with x1 + ... + x3 <= 1 is generated
  // by producing a point with x1 + ... + x3 = 1 and scaling
  // by s, s is distributed according to f(x) = 3x^2 (cdf is F(x)=x^3).
  GMM_random_stochastic_vector3(v);

  s = uniform_real(0,1);
  s = s*s*s;

  for (long i=0; i < 3; i++) {
    v[i] = s*v[i];
  }
  v[3] = 1 - s; 
}


// Produces a random GMM root distribution.
void GMM_random_root(Root &r) {
  GMM_random_stochastic_vector4(r);
}

// Produces a random GMM transition matrix
void GMM_random_edge(TMatrix &tm) {
  GMM_random_stochastic_vector4(tm[0]);
  GMM_random_stochastic_vector4(tm[1]);
  GMM_random_stochastic_vector4(tm[2]);
  GMM_random_stochastic_vector4(tm[3]);
}


void GMM_random_rate_matrix(double tr, TMatrix &Q) {
  std::vector<double> v;
  v.resize(Q.size());

  // Fill in the diagonal
  GMM_random_stochastic_vector4(v);
  for(long i=0; i < 4; i++) {
    Q[i][i] = tr*v[i];
  }

  // Complete the rows.
  long k;
  for(long i=0; i < 4; i++) {
    GMM_random_stochastic_vector3(v);
    k=0;
    for(long j=0; j < 4; j++) {
      if (j == i) continue;
      Q[i][j] = -Q[i][i] * v[k];
      k++;
    }
  }
}

// translating transition matrix to 'our' format
void GMM_matrix_exponential(TMatrix &Q, TMatrix &A) {
  using namespace boost::numeric;

  ublas::matrix<double> QQ(4,4);
  ublas::matrix<double> AA(4,4);

  for(unsigned long i = 0; i < Q.size(); i++) {
    for(unsigned long j = 0; j < Q.size(); j++) {
      QQ(i,j) = Q[i][j];
    }
  }

  AA = expm_pad(QQ);

  for(unsigned long i = 0; i < Q.size(); i++) {
    for(unsigned long j = 0; j < Q.size(); j++) {
      A[i][j] = AA(i,j);
    }
  }
}


// Random GMM transition matrix of a given branch length
void GMM_random_edge_length(double len, TMatrix &tm) {
  double t;
  TMatrix Q, A, B;
  Q.resize(tm.size());
  A.resize(tm.size());
  B.resize(tm.size());
  for(unsigned long l=0; l < tm.size(); l++) {
    Q[l].resize(tm.size());
    A[l].resize(tm.size());
    B[l].resize(tm.size());
  }

  t = uniform_real(-4*len, 0);
  GMM_random_rate_matrix(t, Q);
  GMM_matrix_exponential(Q, A);  

  SSM_random_edge_length(len + t/4, B);


  for(long i = 0; i < 4; i++) {
    for(long j=0; j < 4; j++) {
      tm[i][j] = 0;
      for(long k=0; k < 4; k++) {
        tm[i][j] = tm[i][j] + B[i][k]*A[k][j];
      }
    }
  }
}



// Generates a biologically meaningful GMM matrix (Chang's DLC condition)
void GMM_random_edge_bio_length(double len, TMatrix &tm) {
  TMatrix tmaux;
  long j, k;
  Permutation p;
  p.resize(4);

  tmaux.resize(4);
  for(j=0; j < 4; j++) {
    tmaux[j].resize(4);
  }

  // Loop until there is a permutation that puts it into the DLC form.
  long timeout=0;
  do {
    GMM_random_edge_length(len, tmaux);

    p[0] = max_in_col(tmaux, 0);
    p[1] = max_in_col(tmaux, 1);
    p[2] = max_in_col(tmaux, 2);
    p[3] = max_in_col(tmaux, 3);

    timeout++;
  } while(!(is_permutation(p) && permutation_sign(p) == 1) && (timeout < 1000));
  
  if (timeout >= 1000) {
    std::cout << "ERROR: In sampling for GMM model. Can't generate DLC matrix of length " << len;
    std::cout << std::endl;
    exit(1);
  }

  for(j=0; j < 4; j++) {
    for(k=0; k < 4; k++) {
      tm[j][k] = tmaux[p[j]][k];
    }
  }
}
