/*
 *  sampling.cpp
 *  
 *  Created by Ania M. Kedzierska on 11/11/11.
 *  Copyright 2011 Politecnic University of Catalonia, Center for Genomic Regulation.  This is program can be redistributed, modified or else as given by the terms of the GNU General Public License. 
 *  
 */

#include <iostream>
#include <cstdlib>

#include "random.h"
#include "sampling.h"
#include "miscelania.h"
#include "em.h"

#include "boost/math/distributions/chi_squared.hpp"

// Produces random parameters without length restriction.
void random_parameters(Model &Mod, Parameters &Par){
  long k;
  if (Par.nalpha != Mod.nalpha) {
    std::cout << "Number of states in model does not match parameters";
    exit(1);
  }
  for(k=0; k < Par.nedges; k++) {
    Mod.random_edge(Par.tm[k]);
  }
  Mod.random_root(Par.r);
}


// Produces random parameters of given lengths (as given inside Par).
// The lengths are given as a vector matching the list of edges.
// The generated parameters are biologically meaningful ( DLC as in Chang, 1966), i.e. every diagonal entry is the largest in its column.
void random_parameters_length(Tree &T, Model &Mod, Parameters &Par){
  long k;

  if (Par.nalpha != Mod.nalpha) {
    std::cout << "Number of states in model does not match parameters";
    exit(1);
  }

  if (T.nedges != Par.nedges) {
    std::cout << "Number of edges in tree does not match parameters";
    exit(1);
  }

  for(k=0; k < Par.nedges; k++) {
    Mod.random_edge_bio_length((T.edges[k]).br, Par.tm[k]);
  }
  Mod.random_root(Par.r);
}



// Simulates the data used for testing. Stores the data in data and the parameters in Parsim.
void random_fake_counts(Tree &T, double N, Counts &data, Parameters &Parsim) {

  long i;

  // Fills in simulated counts
  data.N = 0;
  data.c.resize(T.nstleaves);         // assigns memory
  for (i=0; i < T.nstleaves; i++) {   
    data.c[i] = N*joint_prob_leaves(T, Parsim, i);
    data.N = data.N + data.c[i];
  }

  // Copies species names.
  data.species.resize(T.nleaves);     // assigns memory
  for (i=0; i < T.nleaves; i++) {     
    data.species[i] = T.names[i];
  }

  // Fills in the rest of the struct.
  data.nalpha = 4;
  data.nspecies = T.nleaves;
  data.nstates = T.nstleaves;
}


// Puts a random state on st, following the given model and parameters.
void random_state(Tree &T, Model &Mod, Parameters &Par, State &st) {
  long i, e;
  bool updated;
  if (st.len != T.nnodes) {
    std::cout << "Tree and state are not compatible." << std::endl;
  }

  // reset the state to a negative value.
  for(i=0; i < st.len; i++) st.s[i] = -1;

  // random state on the root.
  st.s[T.nleaves] = discrete(Par.r);

  updated = true;
  while(updated) {
    updated = false;
    // run over the edges, and update any states necessary
    for(e=0; e < T.nedges; e++) {
      // if there is something to be updated

      if (st.s[T.edges[e].s] >= 0 && st.s[T.edges[e].t] < 0) {
        updated = true;
        st.s[T.edges[e].t] = discrete(Par.tm[e][st.s[T.edges[e].s]]);
      }
    }
  }
}


// Simulates counts for the full observable and hidden state matrix F.
void random_full_counts(Tree &T, Model &Mod, Parameters &Par, long length, Matrix &F) {
  long i, j, k;
  long a, b;

  State st, sth, stl;

  create_state(st, T.nnodes, T.nalpha);
  create_state(sth, T.nhidden, T.nalpha);
  create_state(stl, T.nleaves, T.nalpha);

  for(i=0; i < T.nstleaves; i++) {
    for(j=0; j < T.nsthidden; j++) {
      F.m[i][j] = 0;
    }
  }
  for (k=0; k < length; k++) {
    random_state(T, Mod, Par, st);

    // the first nleaves nodes are the leaves
    for (i=0; i < T.nleaves; i++) {
      stl.s[i] = st.s[i];
    }

    for (i=T.nleaves; i < T.nnodes; i++) {
      sth.s[i - T.nleaves] = st.s[i];
    }
    a = state2index(stl);
    b = state2index(sth);
    F.m[a][b] = F.m[a][b] + 1.;
  }
}



// Simulates an alignment for a given tree, model and parameters.
void random_data(Tree &T, Model &Mod, Parameters &Par, long length, Alignment &align) {

  if (T.nedges != Par.nedges) {
    std::cout << "Error: Tree and Parameters are not compatible" << std::endl;
  }
  long i, j;
  State st;
  create_state(st, T.nnodes, T.nalpha);

  // The 4  DNA bases:
  align.nalpha = 4;
  align.alpha.resize(align.nalpha);
  align.alpha[0] = 'A';
  align.alpha[1] = 'C';
  align.alpha[2] = 'G';
  align.alpha[3] = 'T';

  align.nspecies = T.nleaves;
  align.len = length;
  align.seq.resize(align.nspecies);
  align.name.resize(align.nspecies);
  for (i=0; i < align.nspecies; i++) {
    (align.seq[i]).resize(length);
    align.name[i] = T.names[i];
  }

  for (i=0; i < length; i++) {
    random_state(T, Mod, Par, st);
    for (j=0; j < align.nspecies; j++) {
      // the first nleaves nodes are the leaves
      align.seq[j][i] = st.s[j];    
    }
  }
}




// Checks if a matrix tm belongs to a given model( for matrices generated during sampling).
bool check_model_matrix(Model &Mod, TMatrix &tm) {
  long i, j, k, l;
  long a, b;
  for(i=0; i < Mod.nalpha; i++) {
    for(j=0; j < Mod.nalpha; j++) {
      for(k=0; k < Mod.nalpha; k++) {
        for(l=0; l < Mod.nalpha; l++) {
          a = Mod.matrix_structure(i, j);
          b = Mod.matrix_structure(k, l);
          if (a == b && tm[i][j] != tm[k][l]) {
            return false;
          }
	}
      }
    }
  }
  return true;
}


bool check_DLC_matrix(TMatrix &tm) {
  long i0;
  for (unsigned long i=0; i < tm[0].size(); i++) {
    i0 = max_in_col(tm, i);
    if (i0 != (long)i) return false;
  }
  return true;
}


bool check_matrix_length(double len, TMatrix &tm) {
  double clen;
  double eps = 1e-10;
  clen = branch_length(tm, tm.size());
  if (boost::math::isnan(clen)) return false;

  if (fabs(clen - len) > eps) return false;
  else return true;
}


bool check_matrix_stochastic(TMatrix &tm) {
  double sum;
  double eps = 1e-10;
  for(unsigned long i=0; i < tm.size(); i++) {
    sum = 0;
    for(unsigned long j=0; j < tm[i].size(); j++) {
      if (tm[i][j] < 0) return false;
      sum = sum + tm[i][j];
    }

    if (fabs(sum - 1) > eps) return false;
  }
  return true;
}
