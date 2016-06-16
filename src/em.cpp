/*
 *  em.cpp
 *  
 *  Created by Ania M. Kedzierska on 11/11/11.
 *  Copyright 2011 Politecnic University of Catalonia, Center for Genomic Regulation.  This is program can be redistributed, modified or else as given by the terms of the GNU General Public License. 
 *  
 */


#include <cstdlib>
#include <iostream>
#include <cmath>

#include "em.h"
#include "matrix.h"
#include "miscelania.h"

////////////////////////////////////////////////////////////////////////////////////////
//
// This block contains joint_prob that compute joint probabilities from the parameters, 
// together with some other functions that use joint_prob.
// This version does not need a precomputed list of states, so it is easier to use,
//  but because it uses index2state very intensively, it is slower than the version below.
//
////////////////////////////////////////////////////////////////////////////////////////


// This computes the joint probability for a choice of states on the hidden nodes (idhid) and
// on the leaves (idleaf).
double joint_prob(Tree &T, Parameters &Par, long idleaf, long idhid) {
  double p;
  long i;
  long a,b;

  // Computes the states and stores them in the auxiliar members of the tree T
  index2state(idhid, T.shidden);
  index2state(idleaf, T.sleaves);

  // root distribution. The root is the first hidden node 
  p = Par.r[T.shidden.s[0]];

  // product of transition matrices
  for(i=0; i < T.nedges; i++) {
    a = T.shidden.s[T.edges[i].s - T.nleaves];      // state on the source node
    if (T.edges[i].t < T.nleaves){                  // if the edge connects to a leaf
      b = T.sleaves.s[T.edges[i].t];                // state on the target leaf
    } else {
      b = T.shidden.s[T.edges[i].t - T.nleaves];    // state on the target node
    }
    p = p*Par.tm[i][a][b];
  }
  return p;
}


// Computes the joint probability on the leaves only, marginalizing for the hidden nodes.
double joint_prob_leaves(Tree &T, Parameters &Par, long idleaf) {
  long i;
  double p;
  p = 0;
  for(i=0; i < T.nsthidden; i++) {
    p = p + joint_prob(T, Par, idleaf, i);
  }
  return p;
}

// Computation of the log-likelihood of the parameters Par given some data on the leaves.
double log_likelihood(Tree &T, Parameters &Par, Counts &data) {
  double l;
  long i;
  if (data.nstates != T.nstleaves) {
    std::cout << "Error: Wrong data length in log_likelihood" << std::endl;
    exit(1);
  }

  l = 0;
  for(i=0; i < T.nstleaves; i++) {     // run over all possible states on the leaves.
    l = l + log(joint_prob_leaves(T, Par, i))*data.c[i];
  }
  return l;
}


// Computes the KL divergence between two sets of parameters. 
double KL_divergence(Tree &T, Parameters &Par1, Parameters &Par2) {
  long i;
  double KL, p1, p2;
  KL = 0;
  for (i=0; i < T.nstleaves; i++) {
    p1 = joint_prob_leaves(T, Par1, i);
    p2 = joint_prob_leaves(T, Par2, i);
    if (p1 == 0 && p2 == 0) continue;      // 0log0 = 0 !!!
    KL = KL + p1*log(p1/p2);
  }
  return KL;
}

// Computes the entropy of a set of parameters.
double entropy(Tree &T, Parameters &Par) {
  long i;
  double E, p;
  E = 0;
  for (i=0; i < T.nstleaves; i++) {
    p = joint_prob_leaves(T, Par, i);
    if (p == 0) continue;      // 0log0 = 0 
    E = E - p*log(p);
  }
  return E;
}

////////////////////////////////////////////////////////////////////////////////////////
//
// Here is a faster version of the functions above, that use a precomputed list of states.
// These ones are the ones used in the critical parts of the EM-algorithm.
// (los estados estan precalculados = rellenar un vector v de forma que v[i][j] es el valor en el nodo j del estado codificado con el número i
//i es un numero de 0 a 4^n-1, y j de 0 a #edges-1_
////////////////////////////////////////////////////////////////////////////////////////


// Faster joint probability with precomputed list of states
double joint_prob_fast(Tree &T, Parameters &Par, StateList &sl, long idleaf, long idhid) {
  double p;
  long i;
  int a,b;

  // root distribution. The root is the first hidden node !
  p = Par.r[sl.h[idhid][0]];

  // product of transition matrices
  for(i=0; i < T.nedges; i++) {
    a = sl.h[idhid][T.edges[i].s - T.nleaves];      // state on the source node
    if (T.edges[i].t < T.nleaves){                  // if the edge connects to a leaf
      b = sl.l[idleaf][T.edges[i].t];               // state on the target leaf
    } else {
      b = sl.h[idhid][T.edges[i].t - T.nleaves];    // state on the target node
    }
    p = p*Par.tm[i][a][b];
  }
  return p;
}


// Fast marginalization of the joint probability to the leaves, using a precomputed
// list of states
double joint_prob_leaves_fast(Tree &T, Parameters &Par, StateList &sl, long idleaf) {
  long i;
  double p;
  p = 0;
  for(i=0; i < T.nsthidden; i++) {
    p = p + joint_prob_fast(T, Par, sl, idleaf, i);
  }
  return p;
}

// Fast log likelihood using a precomputed list of states.
double log_likelihood_fast(Tree &T, Parameters &Par, Counts &data, StateList &sl) {
  double l;
  long i;
  if (data.nstates != T.nstleaves) {
    std::cout << "Error: Wrong data length in log_likelihood" << std::endl;
    exit(1);
  }

  l = 0;
  for(i=0; i < T.nstleaves; i++) {
    l = l + log(joint_prob_leaves_fast(T, Par, sl, i))*data.c[i];
  }
  return l;
}

// Fast KL-divergence using a precomputed list of states.
double KL_divergence_fast(Tree &T, Parameters &Par1, Parameters &Par2, StateList &sl) {
  long i;
  double KL, p1, p2;
  KL = 0;
  for (i=0; i < T.nstleaves; i++) {
    p1 = joint_prob_leaves_fast(T, Par1, sl, i);
    p2 = joint_prob_leaves_fast(T, Par2, sl, i);
    if (p1 == 0 && p2 == 0) continue;      // 0log0 = 0 
    KL = KL + p1*log(p1/p2);
  }
  return KL;
}

// Fast entropy using a precomputed list of states.
double entropy_fast(Tree &T, Parameters &Par, StateList &sl) {
  long i;
  double E, p;
  E = 0;
  for (i=0; i < T.nstleaves; i++) {
    p = joint_prob_leaves_fast(T, Par, sl, i);
    if (p == 0) continue;      // 0log0 = 0 !!!
    E = E - p*log(p);
  }
  return E;
}



// Computes the marginalization to the nodes a, b of the joint probabilities in F.
// Stores the result in the matrix N (nalpha x nalpha). the state of a indexes the rows
// in N, and the state of b indexes the columns.
// The rows of F are indexed by leaf states, and the col's by hidden node states.
// Uses a precomputed list of states for increased speed efficiency.
void two_node_marginalization(Tree &T, Matrix &F, long a, long b, StateList &sl, TMatrix &N) {
  long i,j;
  int sta, stb;

  // initializing N to zero
  for (i=0; i < T.nalpha; i++) {
    for (j=0; j < T.nalpha; j++) {
      N[i][j] = 0;
    }
  }

  for (i=0; i < T.nstleaves; i++) {
    for(j=0; j < T.nsthidden; j++) {

      // Gets the state at the node a
      if (a < T.nleaves) sta = sl.l[i][a];      // if a is a leaf
      else sta = sl.h[j][a - T.nleaves];        // if a is a hidden node

      // Gets the state at the node b
      if (b < T.nleaves) stb = sl.l[i][b];      // if b is a leaf
      else stb = sl.h[j][b - T.nleaves];        // if b is a hidden node

      N[sta][stb] = N[sta][stb] + F.m[i][j];
    }
  }
}

// Computes the marginalization to the node a of the joint probabilities in F.
// Stores the result in the vector s (of length nalpha).
// Uses a precomputed list of states for speed efficiency.
void one_node_marginalization(Tree &T, Matrix &F, long a, StateList &sl, Root &s) {
  long i,j;
  int sta;

  // initializing s to zero
  for (i=0; i < T.nalpha; i++) {
    s[i] = 0;
  }

  for (i=0; i < T.nstleaves; i++) {
    for(j=0; j < T.nsthidden; j++) {

      // Gets the state at the node a
      if (a < T.nleaves) sta = sl.l[i][a];      // if a is a leaf
      else sta = sl.h[j][a - T.nleaves];        // if a is a hidden node

      s[sta] = s[sta] + F.m[i][j];
    }
  }
}



// Uses the EM algorithm to compute the maximum likelihood parameters for the given
// data. Uses the value in Par as starting point.
// eps is an error threshold for stopping the iteration.
// returns the value of the likelihood 
double EMalgorithm(Tree &T, Model &Mod, Parameters &Par, Counts &data, double eps, bool silent) {
  double LikelNew, LikelOld;
  long iter;
  long i, j, e, a, b;
  double pleaf, p;

  TMatrix N;
  Root s;

  Matrix F;
  std::vector<double> br;

  StateList sl;

  // creates a state list for the tree T, and fills it in. This will speed up
  // the most critical loops below.
  create_state_list(sl, T);

  // initializes the big matrix F
  F = create_matrix(T.nstleaves, T.nsthidden);

  // Initializes the auxiliar N and s
  s.resize(T.nalpha);
  N.resize(T.nalpha);
  for (i=0; i < T.nalpha; i++) {
    N[i].resize(T.nalpha);
  }

  // Initial value for likelihood.
  LikelNew = log_likelihood_fast(T, Par, data, sl);
  LikelNew = 0;
  LikelOld = LikelNew + 100.0; 
  iter = 0;

  // The main loop.
  while (fabs(LikelNew - LikelOld) > eps) {
    LikelOld = LikelNew;

    // Performs an E-step. Uses Par to distribute data into the full matrix F.
    for(i=0; i < F.nrows; i++) {

      pleaf = 0;
      for(j=0; j < F.ncols; j++) {
        p = joint_prob_fast(T, Par, sl, i, j);
        pleaf = pleaf + p;
        F.m[i][j] = data.c[i] * p; 
      }

      // now divide by the probability on the leaf.
      for(j=0; j < F.ncols; j++) {
        F.m[i][j] = F.m[i][j] / pleaf;
      }
    }

    // Performs an M-step. Estimates the MLE for the parameters given
    // the full data in F. This is the only model-dependent part.

    // estimate the transition matrices on the edges
    for (e=0; e < T.nedges; e++) {
      a = T.edges[e].s;
      b = T.edges[e].t;

      two_node_marginalization(T, F, a, b, sl, N);
      Mod.mle_edge(N, Par.tm[e]);
    }

    // estimate the root distribution. T.nleaves represents the root node.
    one_node_marginalization(T, F, T.nleaves, sl, s);
    Mod.mle_root(s, Par.r);

    // updates the likelihood
    LikelNew = log_likelihood_fast(T, Par, data, sl);

    iter++;

    if (!silent) {
      std::cout << "  " << iter << ":  L: " << LikelNew;
      std::cout << "  err: " << fabs(LikelNew - LikelOld) << "             \r";

      // forces to print. Otherwise keeps output into a buffer for some time.
      std::cout.flush();   
    }
  }
  if (!silent) {
    std::cout << std::endl;
  }

  delete_matrix(F);
  return LikelNew;
}


// Computes the MLE for the parameters, using counts without hidden nodes.
// The counts are stored in F. Rows are indexed by the leaf states, cols by the hidden states.
// Uses an auxiliar statelist sl to speed up the computation.
void MLE_all_obs(Tree &T, Model &Mod, Parameters &Par, Matrix &F, StateList &sl) {
  long e, i;

  TMatrix N;
  Root s;
  s.resize(T.nalpha);
  N.resize(T.nalpha);
  for (i=0; i < T.nalpha; i++) {
    N[i].resize(T.nalpha);
  }

  for(e=0; e < T.nedges; e++) {
    two_node_marginalization(T, F, T.edges[e].s, T.edges[e].t, sl, N);
    Mod.mle_edge(N, Par.tm[e]);
  }


  one_node_marginalization(T, F, T.nleaves, sl, s);
  Mod.mle_root(s, Par.r);
}



// Sets a fixed initial parameters of JC type.
void initial_parameters(Model &Mod, Parameters &Par) {
  long i,j,k;
  for (k=0; k < Par.nedges; k++) {
    for(i=0; i < Par.nalpha; i++) {
      for(j=0; j < Par.nalpha; j++) {
        if (i==j) Par.tm[k][i][j] = 0.7;
        else Par.tm[k][i][j] = 0.3 / (double)(Par.nalpha-1);
      }
    }
  }
  for (i=0; i < Par.nalpha; i++) {
    Par.r[i] = 1/(double)Par.nalpha;
  }
}
