
/*
 *  model_ssm.cpp
 *  
 *  Created by Ania M. Kedzierska on 11/11/11.
 *  Copyright 2011 Politecnic University of Catalonia, Center for Genomic Regulation.  This is program can be redistributed, modified or else as given by the terms of the GNU General Public License. 
 *  
 */

#include "model_ssm.h"
#include "random.h"
#include <cstdlib>
#include <iostream>
#include <cmath>
#include <algorithm>
#include "miscelania.h"


//////////////////////////////////////////////////////////////////
// Functions specific to the Strand Symmetric model, 4 states
//////////////////////////////////////////////////////////////////


void SSM_matrix(double a, double b, double c, double d, double e, double f, double g, double h, TMatrix &tm) {
  tm[0][0] = a;
  tm[0][1] = b;
  tm[0][2] = c;
  tm[0][3] = d;

  tm[1][0] = e;
  tm[1][1] = f;
  tm[1][2] = g;
  tm[1][3] = h;

  tm[2][0] = h;
  tm[2][1] = g;
  tm[2][2] = f;
  tm[2][3] = e;

  tm[3][0] = d;
  tm[3][1] = c;
  tm[3][2] = b;
  tm[3][3] = a;
}




long SSM_matrix_structure(long i, long j) {
  if (i > 1) {
    i = 3 - i;
    j = 3 - j;
  }  

  if (i == 0 && j == 0) return 0;
  else if (i == 0 && j == 1) return 1;
  else if (i == 0 && j == 2) return 2;
  else if (i == 0 && j == 3) return 3;
  else if (i == 1 && j == 0) return 4;
  else if (i == 1 && j == 1) return 5;
  else if (i == 1 && j == 2) return 6;
  else if (i == 1 && j == 3) return 7;
  else {
    std::cout << "ERROR: unknown condition in K81_matrix_structure." << std::endl;
    exit(1);
  }
}


long SSM_root_structure(long i) {
  if (i == 0 || i == 3) return 0;
  else return 1;
}


Permutation SSM_perm(long a, long b, long c, long d) {
  Permutation p;
  p.resize(4);
  p[0] = a; p[1] = b; p[2] = c; p[3] = d;
  return p;
}

void SSM_list_permutations(std::list<Permutation> &L) {
  L.clear();
  L.push_back(SSM_perm(0,1,2,3));  // det 1
  L.push_back(SSM_perm(0,2,1,3));  // det -1

  L.push_back(SSM_perm(1,0,3,2));  // det 1
  L.push_back(SSM_perm(1,3,0,2));  // det -1

  L.push_back(SSM_perm(2,3,0,1));  // det 1
  L.push_back(SSM_perm(2,0,3,1));  // det -1

  L.push_back(SSM_perm(3,2,1,0));  // det 1
  L.push_back(SSM_perm(3,1,2,0));  // det -1

}




// Computes the mle for the transition matrix of a given edge from the
// the marginalization of the counts on that edge.
// N is a matrix representing the counts on that edge
// tm is the matrix where the mle is stored. tm must be a stochastic matrix
void SSM_mle_edge(TMatrix &N, TMatrix &tm) {
    double s[4];
    long i,j;

    for (i=0; i < 4; i++) {
      s[i] = 0;
      for(j=0; j < 4; j++) {
        s[i] = s[i] + N[i][j];
      }
    }

    tm[0][0] = (N[0][0] + N[3][3])/(s[0]+s[3]);
    tm[0][1] = (N[0][1] + N[3][2])/(s[0]+s[3]);
    tm[0][2] = (N[0][2] + N[3][1])/(s[0]+s[3]);
    tm[0][3] = (N[0][3] + N[3][0])/(s[0]+s[3]);

    tm[1][0] = (N[1][0] + N[2][3])/(s[1]+s[2]);
    tm[1][1] = (N[1][1] + N[2][2])/(s[1]+s[2]);
    tm[1][2] = (N[1][2] + N[2][1])/(s[1]+s[2]);
    tm[1][3] = (N[1][3] + N[2][0])/(s[1]+s[2]);

    tm[2][0] = tm[1][3];
    tm[2][1] = tm[1][2];
    tm[2][2] = tm[1][1];
    tm[2][3] = tm[1][0];

    tm[3][0] = tm[0][3];
    tm[3][1] = tm[0][2];
    tm[3][2] = tm[0][1];
    tm[3][3] = tm[0][0];
}

// Computes the mle for the root distribution from the
// marginalization of the counts to the root node.
// s is a vector representing the counts on the root
// r is a vector where the root distribution is stored.
void SSM_mle_root(Root &s, Root &r) {
  r[0] = 0.5*(s[0] + s[3])/(s[0]+s[1]+s[2]+s[3]);
  r[1] = 0.5*(s[1] + s[2])/(s[0]+s[1]+s[2]+s[3]);
  r[2] = r[1];
  r[3] = r[0];
}


// Produces a random stochastic row.
void SSM_random_row(std::vector<double> &v) {
  double a,b,c,d;

  do {
    a = uniform_real(0, 1);
    b = uniform_real(0, 1);
    c = uniform_real(0, 1);
    d = 1 - a - b - c;
  } while(d < 0);

  v[0] = a;
  v[1] = b;
  v[2] = c;
  v[3] = d;
}


// Computes a random SSM root distribution
void SSM_random_root(Root &r) {
  double x;
  x = uniform_real(0,1);
  r[0] = 0.5*x;
  r[1] = 0.5*(1-x);
  r[2] = 0.5*(1-x);
  r[3] = 0.5*x;
}

// Computes a random ssm matrix with diagonal entries >= 0.7.
void SSM_random_edge(TMatrix &tm) {
  SSM_random_row(tm[0]);
  SSM_random_row(tm[1]);
  
  tm[3][0] = tm[0][3];
  tm[3][1] = tm[0][2];
  tm[3][2] = tm[0][1];
  tm[3][3] = tm[0][0];

  tm[2][0] = tm[1][3];
  tm[2][1] = tm[1][2];
  tm[2][2] = tm[1][1];
  tm[2][3] = tm[1][0];
}


// Random SSM matrix of given length
void SSM_random_edge_length(double len, TMatrix &tm) {
  double a,b,c,d,e,f,g,h;
  double K, v0, lambda, mu, alpha, beta, alphap, betap, s, t, rk, E, F;

  K = exp(-4*len);

  v0 = 1./3.*pow(3.*sqrt(81.*K*K + 3.) + 27.*K, 1./3.) 
     - 1./3.*pow(3.*sqrt(81.*K*K + 3.) - 27.*K, 1./3.);
  s = uniform_real(v0+1, 2);
  rk = (s-1)*(s-1)*(s-1) + (s-1) - 2*K;
  t = uniform_real(0, std::min<double>(2-s, sqrt(rk/(s-1))));

  double gamma;
  gamma = uniform_real(0,1);
  if (gamma > 0.5) t = -t;

  lambda = 0.5*(s+t);
  mu = 0.5*(s-t);

  E = K/(lambda + mu - 1) - (1-lambda)*(1-mu);
  F = K/(lambda + mu - 1) + (1-lambda)*(1-mu);
  alpha = uniform_real(std::max<double>(0, E/lambda), mu);
  beta = uniform_real(std::max<double>(-lambda, E/alpha), std::min<double>(lambda, F/alpha));
  betap = uniform_real(fabs(alpha*beta - K/(lambda + mu - 1))/(1-mu), 1 - lambda);

  gamma = uniform_real(0,1);
  if (gamma > 0.5) betap = -betap;
  
  alphap = (alpha*beta - K/(lambda + mu - 1)) / betap;
  a = 0.5*(lambda + beta);
  b = 0.5*(1 - lambda - betap);
  c = 0.5*(1 - lambda + betap);
  d = 0.5*(lambda - beta);
  e = 0.5*(1 - mu - alphap);
  f = 0.5*(mu + alpha);
  g = 0.5*(mu - alpha);
  h = 0.5*(1 - mu + alphap);


  // Permutation to get rid of f > g
  gamma = uniform_real(0,1);
  if (gamma > 0.5) {
    swap(b,c);   swap(a,d);
    swap(f,g);   swap(e,h);
  }

  // Permutation to get rid of mu + lambda > 1
  gamma = uniform_real(0,1);
  if (gamma > 0.5) {
    swap(a,c);   swap(b,d);
    swap(e,g);   swap(f,h);
  }

  SSM_matrix(a,b,c,d,e,f,g,h, tm);
}



// Random biologically meaningful SSM matrix of given length.
// Biologically meaningful (Chang's DLC condition, diagonal entries are largest in a given column).
void SSM_random_edge_bio_length(double len, TMatrix &tm) {
  TMatrix tmaux;
  long i0, i1, i;

  tmaux.resize(4);
  for(i=0; i < 4; i++) {
    tmaux[i].resize(4);
  }

  // Loop until the permutation that puts into DLC does not change det.
  long timeout = 0;
  do {
    SSM_random_edge_length(len, tmaux);
    i0 = max_in_col(tmaux, 0);
    i1 = max_in_col(tmaux, 1);
    timeout++;
  } while (!((i0 == 0 && i1 == 1) || (i0 == 1 && i1 == 0) || 
             (i0 == 2 && i1 == 3) || (i0 == 3 && i1 == 2)) &&
             (timeout < 1000));

  if (timeout >= 1000) {
    std::cout << "ERROR: In sampling for SSM model. Can't generate DLC matrix of length " << len;
    std::cout << std::endl;
    exit(1);
  }

  SSM_matrix(tmaux[i0][0], tmaux[i0][1], tmaux[i0][2], tmaux[i0][3],
             tmaux[i1][0], tmaux[i1][1], tmaux[i1][2], tmaux[i1][3], tm);  
}
