/*
 *  fisher.cpp
 *  
 *  Created by Ania M. Kedzierska on 11/11/11.
 *  Copyright 2011 Politecnic University of Catalonia, Center for Genomic Regulation.  This is program can be redistributed, modified or else as given by the terms of the GNU General Public License. 
 *  
 */

#include <cstdlib>
#include <cmath>

#include "fisher.h"
#include "matrix.h"
#include "state.h"
#include "state_list.h"
#include "alignment.h"

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/triangular.hpp>
#include <boost/numeric/ublas/lu.hpp>
#include <boost/numeric/ublas/io.hpp>

namespace ublas = boost::numeric::ublas;


//Please refer to the Empar paper for the description of the procedure.

//  The free parameters (without stochastic condition) are encoded as a linear index.
// There will be np*nedges + rdf+1, where
// np: free parameters of the transition matrices
// rdf: df of the root.
// We put the parameters in a matric, first the transition parameters, then the root parameters at the end.
// Root parameters are present only for SSM and GMM.

// Converts an edge-row-column to a linear index identifying the parameter for the given model.
long erc2param(Tree &T, Model &Mod, long e, long r, long c) {
  long p = Mod.matrix_structure(r, c);
  return e*Mod.np + p;
}

// Converts the root distribution index into a parameter index for the model.
long root2param(Tree &T, Model &Mod, long r) {
  return T.nedges*Mod.np + Mod.root_structure(r);
}

// STEP 2
// We then encode the free parameters taking into account the stochastic condition.
// i.e. we remove the diagonal entry and first parameter of the root.


// Converts an edge-row-column to a linear index identifying the free parameter.
long erc2freeparam(Tree &T, Model &Mod, long e, long r, long c) {
  long p = Mod.matrix_structure(r, c);
  long a,b;

  // find first appearence of parameter p.
  for (a=0; a < T.nalpha; a++) {
    for(b=0; b < T.nalpha; b++) {
      if (Mod.matrix_structure(a,b) == p) {
        if (b > a) return e*Mod.df + p - a - 1;
        else if (b < a) return e*Mod.df + p - a;
        else return -1;
        
      }
    }
  }

  return -1;
}

// Converts the root distribution index into a free parameter index.
long root2freeparam(Tree &T, Model &Mod, long r) {
  long p = Mod.root_structure(r);

  if (p == 0) return -1;   // p = 0 is the parameter killed.
  else return T.nedges*Mod.df + p - 1;
}


// Computes the state at the endpoints of an edge.
// T: tree
// e: edge
// sth: state on hidden nodes
// stl: state on leaves

// Output:
// a: the state of the source.
// b: the state of the target.

void edgestate(Tree &T, long e, State &sth, State &stl, long &a, long &b) {
  a = sth.s[T.edges[e].s - T.nleaves];      // state on the source node
  if (T.edges[e].t < T.nleaves){            // if the edge connects to a leaf
    b = stl.s[T.edges[e].t];                // state on the target leaf
  } else {
    b = sth.s[T.edges[e].t - T.nleaves];    // state on the target node
  }
}


void edgestate(Tree &T, long e, std::vector<int> &sth, std::vector<int> &stl, long &a, long &b) {
  a = sth[T.edges[e].s - T.nleaves];      // state on the source node
  if (T.edges[e].t < T.nleaves){          // if the edge connects to a leaf
    b = stl[T.edges[e].t];                // state on the target leaf
  } else {
    b = sth[T.edges[e].t - T.nleaves];    // state on the target node
  }
}


// Returns the root state
long rootstate(Tree &T, State &sth) {
  return sth.s[0];
}

long rootstate(Tree &T, std::vector<int> &sth) {
  return sth[0];
}


// Let k be a parameter index (running over all free parameters on tm + root)
// Every value of k describes a set of states compatible with it (see text).
//   * if k corresponds to a transition matrix on edge e, parameter l, the compatible
//     states are the ones with a state (a,b) on edge e, such that the free parameter
//     on the entry (a,b) in the transition matrix is l.
//   * if k corresponds to a root parameter, then the compatible states are the ones
//     with root state matching that parameter.


// The following functions fill three arrays of probabilities.

// pleaf is a vector of joint probabilities on the leaves
void fill_pleaf(Tree &T, Model &Mod, StateList &sl, Parameters &Par, Array1 &pleaf) {
  long i, j, l, u, v;
  double p;

  // Initializations
  pleaf.resize(T.nstleaves);
  for (i=0; i < T.nstleaves; i++) {
    pleaf[i] = 0;
  }

  // Loop over all states
  for (i=0; i < T.nstleaves; i++) {
    for(j=0; j < T.nsthidden; j++) {

      // Compute leaf probabilities
      p = Par.r[rootstate(T, sl.h[j])];
      for(l=0; l < T.nedges; l++) {
        edgestate(T, l, sl.h[j], sl.l[i], u, v);
        p = p*Par.tm[l][u][v];
      }
      pleaf[i] = pleaf[i] + p;
    }
  }
}



// pcond1 is a matrix of probabilities p(I and k), where "I and k" refers to states with
// leaf state I and compatible with k.
// rows in pcond1 are indexed by leaf states. Columns are indexed by a linear index encoding
// a free parameter.

void fill_pcond1(Tree &T, Model &Mod, StateList &sl, Parameters &Par, Array2 &pcond1) {

  long i,j,l,u,v;
  long k,e,a,b;
  double p;
  long npars = T.nedges*Mod.np + Mod.rdf + 1;

  // Initializations
  pcond1.resize(T.nstleaves);
  for (i=0; i < T.nstleaves; i++) {
    pcond1[i].resize(npars);
    for(k=0; k < npars; k++) {
      pcond1[i][k] = 0;
    }
  }

  // Loop over all states
  for (i=0; i < T.nstleaves; i++) {
    for(j=0; j < T.nsthidden; j++) {
      
      // Loop over edges
      for(e=0; e < T.nedges; e++) {

        // Compute cond1 probabilities
        p = Par.r[rootstate(T, sl.h[j])];
        for(l=0; l < T.nedges; l++) {
          if (l == e) continue;
          edgestate(T, l, sl.h[j], sl.l[i], u, v);
          p = p*Par.tm[l][u][v];
        }

        edgestate(T, e, sl.h[j], sl.l[i], a, b);
        k = erc2param(T, Mod, e, a, b);
        pcond1[i][k] = pcond1[i][k] + p;
      }

      // Take care of root parameters in pcond1
      p = 1;
      for(l=0; l < T.nedges; l++) {
        edgestate(T, l, sl.h[j], sl.l[i], u, v);
        p = p*Par.tm[l][u][v];
      }

      k = root2param(T, Mod, rootstate(T, sl.h[j]));
      pcond1[i][k] = pcond1[i][k] + p;
    }
  }
}



// pcond2 is a 3-dimensional array of probabilities p(I and k and kp) where "I and k and kp" refers
// to states with leaf state I, and compatible with k and kp.

void fill_pcond2(Tree &T, Model &Mod, StateList &sl, Parameters &Par, Array3 &pcond2) {
  long i, j, l, u, v;
  long k, e, a, b;
  long kp, ep, ap, bp;
  double p;

  long npars = T.nedges*Mod.np + Mod.rdf + 1;

  // Initializations
  pcond2.resize(T.nstleaves);
  for (i=0; i < T.nstleaves; i++) {
    pcond2[i].resize(npars);
    for(k=0; k < npars; k++) {
      pcond2[i][k].resize(npars);
      for(l=0; l < npars; l++) {
        pcond2[i][k][l] = 0;
      }
    }
  }

  // Loop over all states
  for (i=0; i < T.nstleaves; i++) {
    for(j=0; j < T.nsthidden; j++) {

      // Loop over pairs of edges for the edge-edge part of pcond2
      for(e=0; e < T.nedges; e++) {
        for(ep=0; ep < T.nedges; ep++) {

          // when e==ep, the corresponding entry must be zero.
          if (e==ep) continue;

          // Compute cond2 probabilities
          p = Par.r[rootstate(T, sl.h[j])];
          for(l=0; l < T.nedges; l++) {
            if (l == e || l == ep) continue;
            edgestate(T, l, sl.h[j], sl.l[i], u, v);
            p = p*Par.tm[l][u][v];
          }

          edgestate(T, e, sl.h[j], sl.l[i], a, b);
          k = erc2param(T, Mod, e, a, b);

          edgestate(T, ep, sl.h[j], sl.l[i], ap, bp);
          kp = erc2param(T, Mod, ep, ap, bp);
      
          pcond2[i][k][kp] = pcond2[i][k][kp] + p;
        }
      }

      // Loop a single edge for the edge-root part of pcond2
      for(e=0; e < T.nedges; e++) {
        // Take care of root parameters in pcond1
        p = 1;
        for(l=0; l < T.nedges; l++) {
          if (l == e) continue;
          edgestate(T, l, sl.h[j], sl.l[i], u, v);
          p = p*Par.tm[l][u][v];
        }

        edgestate(T, e, sl.h[j], sl.l[i], a, b);
        k = erc2param(T, Mod, e, a, b);

        kp = root2param(T, Mod, rootstate(T, sl.h[j]));
  
        pcond2[i][k][kp] = pcond2[i][k][kp] + p;
        pcond2[i][kp][k] = pcond2[i][kp][k] + p;
      }

      // Recall: In pcond2 with two root parameters k and kp,
      // we always have p(I and k and kp) = 0.
    }
  }  
}


// Puts the parameter values in a vector indexed as in the Fisher matrix for the model.
// (4 parameters per edge on K81).
void get_param_vector(Tree &T, Model &Mod, Parameters &Par, std::vector<double> &param) {
  long e, a, b, k;
  long npar = T.nedges * Mod.np + Mod.rdf + 1;
  param.resize(npar);

  for(e=0; e < T.nedges; e++) {
    for(a=0; a < T.nalpha; a++) {
      for(b=0; b < T.nalpha; b++) {
        k = erc2param(T, Mod, e, a, b);
        param[k] = Par.tm[e][a][b];
      }
    }
  }

  for(a=0; a < T.nalpha; a++) {
    k = root2param(T, Mod, a);
    param[k] = Par.r[a];
  }
}

// Gets the free parameters for a branch
void get_branch_free_param_vector(Tree &T, Model &Mod, Parameters &Par, long br, std::vector<double> &param) {
  long a, b, k;
  param.resize(Mod.df);

  for(a=0; a < T.nalpha; a++) {
    for(b=0; b < T.nalpha; b++) {
      k = erc2freeparam(T, Mod, 0, a, b);
      if(k < 0) continue;
      param[k] = Par.tm[br][a][b];
    }
  }
}


// Puts the parameter values in a vector indexed as in the Fisher matrix for the free parameters.
// (3 parameters per edge on K81).
void get_free_param_vector(Tree &T, Model &Mod, Parameters &Par, std::vector<double> &param) {
  long e, a, b, k;
  long npar = T.nedges * Mod.df + Mod.rdf;
  param.resize(npar);

  for(e=0; e < T.nedges; e++) {
    for(a=0; a < T.nalpha; a++) {
      for(b=0; b < T.nalpha; b++) {
        k = erc2freeparam(T, Mod, e, a, b);
        if(k < 0) continue;
        param[k] = Par.tm[e][a][b];
      }
    }
  }

  for(a=0; a < T.nalpha; a++) {
    k = root2freeparam(T, Mod, a);
    if(k < 0) continue;
    param[k] = Par.r[a];
  }
}





// Computes the observed Fisher information of all Parameters in Par, on the model with
// hidden nodes and stores them in the matrix Imod ( no approximation, so the matrix is not diagonal). 
// The rows and columns of Imod are indexed by a number between 0 and np*nedges + rdf + 1, encoding the free parameters
// of the model. The stochastic condition is not taken into account.

// If data is not NULL, computes observed Fisher instead of expected Fisher.

// The formula for the Fisher info of the model parameters (without stochastic constraints) is:
// I(k, k') = sum_{leaf states I} N p(I and k)*p(I and k') / (p_k*p_k'*p(I))
//                            - p(I and k and k') / (p_k*p_k')   <----- This term only when k!=k'
// Here p_k denotes the value of the parameter k.

// The observed Fisher info, is given by the same formula, but inside the sum must add the factor
// x_I / N*p(I). 

void Fisher_information_model(Tree &T, Model &Mod, Parameters &Par, long N, Counts *data, Array2 &Imod) {
  long i, j, k;

  double factor;

  std::vector<double> param;

  Array1 pleaf;
  Array2 pcond1;
  Array3 pcond2;

  StateList sl;

  long npars = Mod.np*T.nedges + Mod.rdf + 1;

  Array2 Ineg;
  Ineg.resize(npars);
  Imod.resize(npars);
  for (i=0; i < npars; i++) {
    Ineg[i].resize(npars);
    Imod[i].resize(npars);
    for(j=0; j < npars; j++) {
      Ineg[i][j] = 0;
      Imod[i][j] = 0;
    }
  }

  if (data != NULL && data->nspecies != T.nleaves) {
    std::cout << "ERROR: In Fisher_information_model. Counts don't match the tree." << std::endl;
    exit(-1);
  }

  create_state_list(sl, T);

  // Fills stuff
  fill_pleaf(T, Mod, sl, Par, pleaf);
  fill_pcond1(T, Mod, sl, Par, pcond1);
  fill_pcond2(T, Mod, sl, Par, pcond2);

  get_param_vector(T, Mod, Par, param);

  // Fills entries of Imod one by one.
  for (i=0; i < npars; i++) {
    for(j=0; j < npars; j++) {
      for(k=0; k < T.nstleaves; k++) {
        if (data == NULL) {    // Expected Fisher
          factor = (double) N;
        } else {               // Observed Fisher
          factor = data->c[k] / pleaf[k];
        }

        Imod[i][j] = Imod[i][j] +
          factor * pcond1[k][i] * pcond1[k][j] / pleaf[k];

        Ineg[i][j] = Ineg[i][j] + factor * pcond2[k][i][j];
      }

      // We add up the negatives apart to minimize numerical errors.
      Imod[i][j] = Imod[i][j] - Ineg[i][j];

    }
  }
}



// Takes The Fisher info matrix corresponding to a model Imod, and outputs a Fisher info matrix
// for the free parameters. This means, taking into account the stochastic condition. What we do
// is replace the diagonal parameter by "1 - sum of other parameters"
// Imod has np*nedges + rdf + 1 rows and columns
// Ifree has df*nedges + rdf rows and columns.

// Here np means number of parameters, while df means degrees of freedom. For
// example, K81 has np = 4 and df = 3. SSM has np = 8 and df = 6. The fact that SSM has two 
// stochastic constraints instead of one, means that SSM must be treated separately.

void Fisher_information_free(Tree &T, Model &Mod, Array2 &Imod, Array2 &Ifree) {
  long i, j, e, p, l;

  long l1, l2;

  long nfreepar = Mod.df*T.nedges + Mod.rdf;
  long nmodpar = Mod.np*T.nedges + Mod.rdf + 1;

  Ifree.resize(nfreepar);
  for(i=0; i < nfreepar; i++) {
    Ifree[i].resize(nfreepar);
    for(j=0; j < nfreepar; j++) {
      Ifree[i][j] = 0;
    }
  }

  // Counts the number of times a parameter appears in the the transition matrices and root dist.
  std::vector<double> coeff;
  coeff.resize(nmodpar);
  for(i=0; i < nmodpar; i++) {
    coeff[i] = 0;
  }

  for(e=0; e < T.nedges; e++) {
    for(i=0; i < T.nalpha; i++) {
      for(j=0; j < T.nalpha; j++) {
        p = erc2param(T, Mod, e, i, j);
        coeff[p] = coeff[p] + 1;  // adds one to the corresponding parameter index.
      }
    }
  }
  for (i=0; i < T.nalpha; i++) {
    p = root2param(T, Mod, i);
    coeff[p] = coeff[p] + 1;
  }

  // translation between free parameters and model parameters.
  std::vector<long> modpar;     // corresponding model parameter index
  std::vector<long> modparkill; // corresponding model parameter which is killed.

  modpar.resize(nfreepar);
  modparkill.resize(nfreepar);

  // Transition matrices
  for(e=0; e < T.nedges; e++) {
    for(i=0; i < T.nalpha; i++) {
      for(j=0; j < T.nalpha; j++) {
        l = erc2freeparam(T, Mod, e, i, j);
        if (l < 0) continue;   // negative l means the parameter is killed.
        modpar[l] = erc2param(T, Mod, e, i, j);
        modparkill[l] = erc2param(T, Mod, e, i, i);
      }
    }
  }

  // Root
  if (Mod.rdf > 0) {
    for(i=0; i < T.nalpha; i++) {
      l = root2freeparam(T, Mod, i);
      if (l < 0) continue;   // negative l means the parameter is killed.
      modpar[l] = root2param(T, Mod, i);
      modparkill[l] = root2param(T, Mod, 0);
    }
  }

  double c1, c2;
  for(l1=0; l1 < nfreepar ; l1++) {
    c1 = coeff[modpar[l1]]/coeff[modparkill[l1]];
    for(l2=0; l2 < nfreepar; l2++) {
      c2 = coeff[modpar[l2]]/coeff[modparkill[l2]];

      Ifree[l1][l2] = Imod[modpar[l1]][modpar[l2]]
           - c1 * Imod[modparkill[l1]][modpar[l2]]
           - c2 * Imod[modpar[l1]][modparkill[l2]]
           + c1*c2 * Imod[modparkill[l1]][modparkill[l2]];
    }
  }
}


// Computes the fisher information for the model with hidden nodes, and stores the result in I.
// I is a matrix with df*nedges + rdf rows and cols.
// the index for the rows and columns encodes a free parameter. 
// 0: first free param on edge 0
// 1: second free param on edge 0
// ...
// df: first free param on edge 1
// ...

void Fisher_information(Tree &T, Model &Mod, Parameters &Par, long N, Array2 &I) {
  Array2 Imod;

  Fisher_information_model(T, Mod, Par, N, NULL, Imod);
  Fisher_information_free(T, Mod, Imod, I);
}


void Observed_Fisher_information(Tree &T, Model &Mod, Parameters &Par, Counts &data, Array2 &I) {
  Array2 Imod;

  Fisher_information_model(T, Mod, Par, data.N, &data, Imod);
  Fisher_information_free(T, Mod, Imod, I);
}



/* Matrix inversion routine.
    Uses lu_factorize and lu_substitute in uBLAS to invert a matrix */
template<class T>
bool matrix_inverse (const ublas::matrix<T>& input, ublas::matrix<T>& inverse) {
  using namespace boost::numeric::ublas;
  typedef permutation_matrix<std::size_t> pmatrix;
  // create a working copy of the input
  matrix<T> A(input);
  // create a permutation matrix for the LU-factorization
  pmatrix pm(A.size1());
  // perform LU-factorization
  int res = lu_factorize(A,pm);
  if( res != 0 ) return false;
  // create identity matrix of "inverse"
  inverse.assign(ublas::identity_matrix<T>(A.size1()));
  // backsubstitute to get the inverse
  lu_substitute(A, pm, inverse);
  return true;
}


// Computes the full covariance matrix for the MLE.
void full_MLE_covariance_matrix(Tree &T, Model &Mod, Parameters &Par, long N, Array2 &Cov) {
  long i, j;

  long npars = Mod.df*T.nedges + Mod.rdf;
  ublas::matrix<double> CC(npars, npars);
  ublas::matrix<double> II(npars, npars);
  Array2 I;

  Fisher_information(T, Mod, Par, N, I);    

  // Need to convert from Matrix to ublas::matrix<double>
  for(i=0; i < npars; i++) {
    for(j=0; j < npars; j++) {
      II(i, j) = I[i][j];
    }
  }

  bool res = matrix_inverse(II, CC);
  if (!res) {
    std::cout << "Could not invert the Fisher information matrix." << std::endl;
    exit(1);
  }

  Cov.resize(npars);
  for(i=0; i < npars; i++) {
    Cov[i].resize(npars);
    for(j=0; j < npars; j++) {
      Cov[i][j] = CC(i,j);
    }
  }
}


// Computes the full covariance matrix for the MLE, using observed Fisher info.
void full_MLE_observed_covariance_matrix(Tree &T, Model &Mod, Parameters &Par, Counts &data, Array2 &Cov) {
  long i, j;

  long npars = Mod.df*T.nedges + Mod.rdf;
  ublas::matrix<double> CC(npars, npars);
  ublas::matrix<double> II(npars, npars);
  Array2 I;

  Observed_Fisher_information(T, Mod, Par, data, I);    

  // Need to convert from Matrix to ublas::matrix<double>
  for(i=0; i < npars; i++) {
    for(j=0; j < npars; j++) {
      II(i, j) = I[i][j];
    }
  }

  bool res = matrix_inverse(II, CC);
  if (!res) {
    std::cout << "Could not invert the Fisher information matrix." << std::endl;
    exit(1);
  }

  Cov.resize(npars);
  for(i=0; i < npars; i++) {
    Cov[i].resize(npars);
    for(j=0; j < npars; j++) {
      Cov[i][j] = CC(i,j);
    }
  }
}



// Extracts the covariance matrix for a given branch from the full covariance matrix.
void branch_covariance_matrix(Model &Mod, Array2 &Covfull, long b, Array2 &Covbr) {
  long i,j;

  Covbr.resize(Mod.df);
  for(i=0; i < Mod.df; i++) {
    Covbr[i].resize(Mod.df);
    for(j=0; j< Mod.df; j++) {
      Covbr[i][j] = Covfull[Mod.df*b + i][Mod.df*b + j];
    }
  }
}


// Extracts the inverted covariance matrix for a given branch from the full covariance matrix.
void branch_inverted_covariance_matrix(Model &Mod, Array2 &Covfull, long b, Array2 &Covbri) {
  long i,j;

  ublas::matrix<double> CC(Mod.df, Mod.df);
  ublas::matrix<double> II(Mod.df, Mod.df);

  for(i=0; i < Mod.df; i++) {
    for(j=0; j< Mod.df; j++) {
      CC(i,j) = Covfull[Mod.df*b + i][Mod.df*b + j];
    }
  }

  bool res = matrix_inverse(CC, II);
  if (!res) {
    std::cout << "Could not invert the Covariance matrix." << std::endl;
    exit(1);
  }

  Covbri.resize(Mod.df);
  for(i=0; i < Mod.df; i++) {
    Covbri[i].resize(Mod.df);
    for(j=0; j< Mod.df; j++) {
      Covbri[i][j] = II(i,j);
    }
  }
}

