/*
 *  parameter_test.cpp
 *  
 *  Created by Ania M. Kedzierska on 11/11/11.
 *  Copyright 2011 Politecnic University of Catalonia, Center for Genomic Regulation.  This is program can be redistributed, modified or else as given by the terms of the GNU General Public License. 
 *  
 */
// as in parameter.cpp, this file contains additional functions I used for testing and getting a feel of the problem.

#include "parameter_tests.h"
#include "em.h"
#include "sampling.h"
#include "parameters.h"
#include "miscelania.h"
#include "random.h"
#include "alignment.h"
#include "fisher.h"

#include <vector>
#include <algorithm>
#include "boost/math/distributions/chi_squared.hpp"
#include "boost/math/special_functions/fpclassify.hpp"
#include "boost/math/special_functions/erf.hpp"
#include <boost/math/special_functions/beta.hpp>

#include <cmath>
#include<iostream>
#include <fstream>
#include <list>


// Those values should be read from a file, here they were calculated prior to the simulations
double get_scale_constant(Model &Mod) {
  if (Mod.m == JC) return 11.8877007085245;
  else if (Mod.m == K80) return 14.3316208140274;
  else if (Mod.m == K81) return 17.3950207294731;
  else if (Mod.m == SSM) return 17.3738070524995;
  else return 0;
}

// multinomial sampling	

// p: vector with probabilities
// N: length of the sampled vector
// x: output vector with a single multinomial sample of length p.size(), 
// the entry x[s] is the count for the state s and these entries add up to N.

void multinomial_sample(std::vector<double> &p, long N, std::vector<double> &x) {
	unsigned long i;
	long s;
	x.resize(p.size(), 0);   // counts for every possible state. Initially set to 0.
	for (i=0; i < p.size(); i++) {
		x[i] = 0;
	}
	
	for(long j=0; j < N; j++) {
		s = discrete(p);
		//std::cout << "sample: " << s << "             \n";
		
		x[s]++;
	}
	
}


// sum of KL's for the first row of the transition matrix.
void KL_divergence_edges(Tree &T, Model &Mod, Parameters &Par1, Parameters &Par2, std::vector<double> &KL) {
  long e, i, r, rr;
  double d, dacc;
  r = 0;
  rr = 1;

	//if unimplemented modelwas called
  if (Mod.m != JC && Mod.m != K80  && Mod.m != K81  && Mod.m != SSM) {
    std::cout << "ERROR: KL_divergence_edges not implemented for this model." << std::endl;
    exit(-1);
  }

  KL.resize(T.nedges);
  for (e=0; e < T.nedges; e++) {
    dacc = 0;
    for (i=0; i < T.nalpha; i++) {
      if (Par1.tm[e][r][i] > 0) {
        d = Par1.tm[e][r][i]*log(Par1.tm[e][r][i]/Par2.tm[e][r][i]);
      } else {
        d = 0;
      }
      dacc = dacc + d;
    } 

    if (Mod.m == SSM) {
      for (i=0; i < T.nalpha; i++) {
        if (Par1.tm[e][rr][i] > 0) {
          d = Par1.tm[e][rr][i]*log(Par1.tm[e][rr][i]/Par2.tm[e][rr][i]);
        } else {
          d = 0;
        }
        dacc = dacc + d;
      }
    }
    KL[e] = dacc;
  }
}


// Computes a chi^2 statistic. Covbr is the inverse covariance matrix for a branch.
double chi2_mult(std::vector<double> &mu, std::vector<double> &x, Array2 &Covbr){
  double X = 0.0;
  for(unsigned int i=0; i < mu.size(); i++) {
    for(unsigned int j=0; j < mu.size(); j++) {		
      X = X + Covbr[j][i]*(x[i] - mu[i])*(x[j] - mu[j]);
    }
  }
  return X;	
}
	
	
bool double_pointer_comparison(const double *a, const double *b) {
  return *a < *b;
}

// version 2:
void BH(std::vector<double> &pvals, std::vector<double> &qvals) {
	long i, n;

        // Copies pvals to qvals
        for(unsigned long j=0; j < pvals.size(); j++) {
	  qvals[j] = pvals[j];
        }

     	std::vector<double*> pvals_p;        // the vector of pointers
	n = (long) pvals.size();

        pvals_p.resize(n);
        for(long j=0; j < n; j++) {
	     pvals_p[j] = &qvals[j];
        }

        // Sorts the vector of pointers only.
	std::sort(pvals_p.begin(), pvals_p.end(), double_pointer_comparison);

	     // Here corrects p-values using the sorted vector of pointers pval_p.
        // *(pvals_p[i]) is the value of the i-th smallest number in qvals,
        // moreover, modifying its value here, changes it in the vector qvals.
	
	for(i=n-2; i >= 0; i--) {
	  *(pvals_p[i]) = *(pvals_p[i]) * ((double)n/(double)(i+1));
	  if (*(pvals_p[i]) > *(pvals_p[i+1])) *(pvals_p[i]) = *(pvals_p[i+1]);
	}
}

// Computes the p-value of X for a chi^2 distribution with a given df
double pvalue_chi2 (double X, int df) {
	if (boost::math::isinf(X)) {
		return 0.;
	} else {
		boost::math::chi_squared chi2(df);
		return 1 - boost::math::cdf(chi2, X);
	}
}

// Uses Fischer's method to combine the p-values
double Fisher_combined_pvalue(std::vector<double> &pvals) {
  int i;
  long n;
  double X2;
  n = pvals.size();

  // We compute the Fisher's statistic
  X2 = 0;
  for (i=0; i < n; i++) { 
    X2 = X2 + -2*log(pvals[i]);
  }

  // Return the p-value of Fisher's statistc
  return pvalue_chi2(X2, 2*n);  
}


// Normalizes a vector of data.
void normalize(std::vector<double> &x) {
  unsigned long i;
  double mean, std;
  mean = 0.;
  for (i = 0; i < x.size(); i++) {
    mean = mean + x[i];
  }
  mean = mean / (double)x.size();

  std = 0.;
  for (i=0; i < x.size(); i++) {
    std = std + (x[i] - mean)*(x[i] - mean);
  }
  std = sqrt(std);

  for (i=0; i < x.size(); i++) {
    x[i] = (x[i] - mean)/std;
  }
}


// Combines p-values using Z-score.
double Zscore_combined_pvalue(std::vector<double> &pvals) {

  unsigned long i;
  double z, zcomb, pcomb;
  
  zcomb = 0.;
  for (i=0; i < pvals.size(); i++) {

    // If some p-value is 0, the score will be 0.
    if (pvals[i] <= 0) {
      return 0; 
    }

    // If some p-value is 1, the score will be 1.
    if (pvals[i] >= 1) {
      return 1; 
    }

    z = sqrt(2)*boost::math::erf_inv(2*(1-pvals[i]) - 1);       // Sure ???
    zcomb = zcomb + z;
  }
  zcomb = zcomb / sqrt((double) pvals.size());
  pcomb = 1 - 0.5*(1+erf(zcomb/sqrt(2)));
  return pcomb;
}



// Performs N repetitions of the EM algorithm and computes combined p-values for every edge.
// Parameters:
// Nrep:    Number of repetitions
// length:  Length of the alignment of simulated data
// pvals:   Vector of p-values indexed by the edges in T
// data_prefix: If different from != "" stores the output data 

void parameter_test(Tree &T, Model &Mod, long Nrep, long length, double eps, std::vector<double> &pvals, std::string data_prefix, bool save_mc_exact){

  long iter;
  long i, r;

  double df, C;
  double distance, KL;
	KL=0;
	distance=0;
  double likel;
   
	
  Parameters Parsim, Par, Par_noperm;
  Alignment align;
  Counts data;

  double eps_pseudo = 0.001;     // Amount added to compute the pseudo-counts.

  StateList sl;

  bool save_data = (data_prefix != "");
  
	
  std::string output_filename;
  std::stringstream output_index;
  std::ofstream logfile;
  std::ofstream logdistfile;

  std::ofstream out_chi2;
  std::ofstream out_br;
  std::ofstream out_brPerc;

  std::ofstream out_pvals;
  std::ofstream out_pvals_noperm;
  std::ofstream out_qvals;
  std::ofstream out_bound;
  std::ofstream out_variances;
  std::ofstream out_qvalsComb;
  std::ofstream out_qvalsCombzscore;
  std::ofstream out_covmatrix;

  std::ofstream out_parest;
  std::ofstream out_parsim;

  std::vector<double> KLe;
  std::vector<std::vector<double> > chi2_array; // an array of chi2 for every edge.
  std::vector<std::vector<double> > mult_array; // an array of mult for every edge.
  std::vector<std::vector<double> > br_array; // an array of br. length for every edge.
  std::vector<std::vector<double> > br_arrayPerc; // an array of br. length for every edge.

  std::vector<std::vector<double> > cota_array; // an array of upper bounds of the diff in lengths for every edge.
  std::vector<std::vector<double> > pval_array; // an array of pvals for every edge.
  std::vector<std::vector<double> > pval_noperm_array;
  std::vector<std::vector<double> > qval_array; // an array of qvalues for every edge.	
  std::vector<std::vector<double> > variances_array; // an array of theoretical variances.
  std::vector<std::vector<double> > parest_array; // array of estimated parameters
  std::vector<std::vector<double> > parsim_array; // array of simulation parameters

	//  ci_binom ci_bin; // condfidence interval
  std::vector<std::vector<ci_binom> > CIbinomial ; //  	vector of CIs
	
  std::list<long> produced_nan;

  long npars = T.nedges*Mod.df + Mod.rdf;

  // Initializing pvals
  pvals.resize(T.nedges);

  // Initialize the parameters for simulation of K81 data for testing
  create_parameters(Par, T);
  create_parameters(Parsim, T);

  // Initializing data structures
  KLe.resize(T.nedges);

	pval_array.resize(T.nedges);
        pval_noperm_array.resize(T.nedges);
        qval_array.resize(T.nedges);
	chi2_array.resize(T.nedges);
	mult_array.resize(T.nedges);
	br_array.resize(T.nedges);
        br_arrayPerc.resize(T.nedges);
	cota_array.resize(T.nedges);
        variances_array.resize(npars);
        parest_array.resize(npars);
        parsim_array.resize(npars);

	// initialize to 0's
  for (i=0; i < T.nedges; i++) {
      pval_array[i].resize(Nrep, 0);
      pval_noperm_array[i].resize(Nrep, 0);
      qval_array[i].resize(Nrep, 0);
          chi2_array[i].resize(Nrep, 0);
	  mult_array[i].resize(Nrep, 0);
	  br_array[i].resize(Nrep, 0);
          br_arrayPerc[i].resize(Nrep, 0);
	  cota_array[i].resize(Nrep, 0);
  }

  for(i=0; i < npars; i++) {
    variances_array[i].resize(Nrep, 0);
    parest_array[i].resize(Nrep, 0);
    parsim_array[i].resize(Nrep, 0);
  }

  // Information about the chi^2.
  df = Mod.df;
  C = get_scale_constant(Mod);


  if (save_data) {
    logfile.open((data_prefix + ".log").c_str(), std::ios::out);
    logfile << "model:  " << Mod.name << std::endl;
    logfile << "length: " << length << std::endl;
    logfile << "eps:    " << eps << std::endl;
    logfile << "nalpha: " << T.nalpha << std::endl;
    logfile << "leaves: " << T.nleaves << std::endl;
    logfile << "tree:   " << T.tree_name << std::endl;
    logfile << std::endl;
    logdistfile.open((data_prefix + ".dist.log").c_str(), std::ios::out);

    out_chi2.open(("out_chi2-" + data_prefix + ".txt").c_str(), std::ios::out);
    out_br.open(("out_br-" + data_prefix + ".txt").c_str(), std::ios::out);
    out_brPerc.open(("out_brPerc-" + data_prefix + ".txt").c_str(), std::ios::out);
    out_pvals.open(("out_pvals-" + data_prefix + ".txt").c_str(), std::ios::out);
    out_pvals_noperm.open(("out_pvals_noperm-" + data_prefix + ".txt").c_str(), std::ios::out);
    out_qvals.open(("out_qvals-" + data_prefix + ".txt").c_str(), std::ios::out);
    out_variances.open(("out_variances-" + data_prefix + ".txt").c_str(), std::ios::out);
    out_parest.open(("out_params-est-" + data_prefix + ".txt").c_str(), std::ios::out);
    out_parsim.open(("out_params-sim-" + data_prefix + ".txt").c_str(), std::ios::out);
    out_bound.open(("out_bound-" + data_prefix + ".txt").c_str(), std::ios::out);
    out_qvalsComb.open(("out_qvalsComb-" + data_prefix + ".txt").c_str(), std::ios::out);
    out_qvalsCombzscore.open(("out_qvalsCombzscore-" + data_prefix + ".txt").c_str(), std::ios::out);

    out_parsim.precision(15);
    out_parest.precision(15);
    out_variances.precision(15);
  }
		
		
	// uncomment the 2 following lines if want to fix the parameters
	// random_parameters_length(T, Mod, Parsim);  
	 //random_data(T, Mod, Parsim, length, align);

  for (iter=0; iter < Nrep; iter++) {
    std::cout << "iteration: " << iter << "             \n";

    // Produces an alignment from random parameters
    random_parameters_length(T, Mod, Parsim);  
	  
	  
    random_data(T, Mod, Parsim, length, align);
    get_counts(align, data);
    add_pseudocounts(eps_pseudo, data);

    // Saving data
    if (save_data) {
      output_index.str("");
      output_index << iter;
      output_filename = data_prefix + "-" + output_index.str(); 
      save_alignment(align, output_filename + ".fa");
      save_parameters(Parsim, output_filename + ".sim.dat");
    }

    // random initial conditions
    //    random_parameters(Mod, Par);

    // fixed initial conditions
    initial_parameters(Mod, Par);

    // Runs the EM
    likel = EMalgorithm(T, Mod, Par, data, eps);

    // If algorithm returns NaN skip this iteration.
    if (boost::math::isnan(likel)) {
      produced_nan.push_back(iter);
      continue;
    }
    copy_parameters(Par, Par_noperm);

    // Chooses the best permutation.
    guess_permutation(T, Mod, Par);
 
    distance = parameters_distance(Parsim, Par);
    // tests using KL score:    
	// KL = KL_divergence_fast(T, Parsim, Par, sl);
    // KL = KL_divergence(T, Par, Parsim);
	
	 
      // estimated counts: Par ; original: Parsim
      std::vector<double> counts_est;
      counts_est.resize(T.nalpha, 0);

		
      // calculate the cov matrix
      std::vector<std::vector<double> > Cov;
      Array2 Cov_br;

      full_MLE_covariance_matrix(T, Mod, Parsim, length, Cov);

      if(save_data) {
          save_matrix(Cov, output_filename + ".cov.dat");   	
      }

      // Save the covariances in an array
      std::vector<double> param;
      std::vector<double> param_sim;

      param.resize(npars);
      param_sim.resize(npars);

      get_free_param_vector(T, Mod, Par, param);
      get_free_param_vector(T, Mod, Parsim, param_sim);

      for(i=0; i < npars; i++) {
	variances_array[i][iter] = Cov[i][i];
        parsim_array[i][iter] = param_sim[i];
        parest_array[i][iter] = param[i];
      }

      std::vector<double> xbranca, xbranca_noperm, mubranca;
      double chi2_noperm;
      xbranca.resize(Mod.df);
      xbranca_noperm.resize(Mod.df);
      mubranca.resize(Mod.df);
      for (i=0; i < T.nedges; i++) {
		  r = 0; // row to be fixed
		  // Extracts the covariance matrix, 1 edge
		  
				  branch_inverted_covariance_matrix(Mod, Cov, i, Cov_br);
                  get_branch_free_param_vector(T, Mod, Parsim, i, mubranca);
                  get_branch_free_param_vector(T, Mod, Par, i, xbranca);
                  get_branch_free_param_vector(T, Mod, Par_noperm, i, xbranca_noperm);
			  
			      chi2_array[i][iter] = chi2_mult(mubranca, xbranca, Cov_br);	
                  chi2_noperm = chi2_mult(mubranca, xbranca_noperm, Cov_br);
		  
		
                  pval_array[i][iter] =  pvalue_chi2(chi2_array[i][iter], Mod.df);  
                  pval_noperm_array[i][iter] = pvalue_chi2(chi2_noperm, Mod.df);             	  

			      br_array[i][iter] = T.edges[i].br - branch_length(Par.tm[i], T.nalpha);
				  br_arrayPerc[i][iter] = branch_length(Par.tm[i], T.nalpha)/T.edges[i].br;


		// Upper bound on the parameter distance using multinomial:
		//  cota_array[i][iter] = bound_mult(Parsim.tm[i], Xm, length); 
	    // and using the  L2 bound
		  cota_array[i][iter] = branch_length_error_bound_mult(Parsim.tm[i], Par.tm[i]);  
		  out_br <<  br_array[i][iter]  << " ";
		  out_brPerc <<  br_arrayPerc[i][iter]  << " ";
		  
		  out_bound  <<  cota_array[i][iter] << " ";
		  out_chi2 << chi2_array[i][iter] << " ";
     			
			}
         out_chi2 << std::endl;		
		 out_bound  <<  std::endl; 
		 out_br << std::endl;
		 out_brPerc << std::endl;
 


	
    // Saves more data.
    if (save_data) {
      logfile << iter << ": " << distance << "   " << KL << std::endl;
      save_parameters(Par, output_filename + ".est.dat");

      logdistfile << iter << ": ";
      logdistfile << parameters_distance_root(Par, Parsim) << " ";
      for(int j=0; j < T.nedges; j++) {
        logdistfile << parameters_distance_edge(Par, Parsim, j) << " ";
      }
      logdistfile << std::endl;
    }
  
} // close iter loop here

  // Correct the p-values
  for(i=0; i < T.nedges; i++) {
       BH(pval_array[i], qval_array[i]);
	//save them
  }

  if (save_mc_exact) {
    for(long iter=0; iter < Nrep; iter++) {
      for(long i=0; i < T.nedges; i++) {
        out_pvals << pval_array[i][iter] << "  ";
        out_pvals_noperm << pval_noperm_array[i][iter] << "  ";
        out_qvals << qval_array[i][iter] << "  ";
      }
      out_pvals  << std::endl;
      out_pvals_noperm << std::endl;
      out_qvals  << std::endl;

      for(long i=0; i < npars; i++) {
        out_variances << variances_array[i][iter] << "  ";
        out_parsim << parsim_array[i][iter] << "  ";
        out_parest << parest_array[i][iter] << "  ";
      }
      out_variances << std::endl;
      out_parsim << std::endl;
      out_parest << std::endl;
    }
  }		
	
	// now combine the pvalues
   for(i=0; i < T.nedges; i++) {		
	pvals[i] = Fisher_combined_pvalue(pval_array[i]);
    //using the Zscore it goes like this: pvals[i] = Zscore_combined_pvalue(pval_array[i]);
	if (save_mc_exact) {	   out_qvalsComb <<  pvals[i] << "  " ;
	out_qvalsCombzscore << Zscore_combined_pvalue(pval_array[i]) << " "; 
	}
  } 

  // Close files
  if (save_data) {
    logdistfile.close();
    logfile.close();
  }

if (save_mc_exact) {
	out_chi2.close();		
	out_bound.close();
        out_variances.close();
        out_parest.close();
        out_parsim.close();
	out_br.close();
	out_brPerc.close();

	out_pvals.close();
        out_qvals.close();
	out_qvalsComb.close();
	out_qvalsCombzscore.close();
        out_covmatrix.close();
	}

  // Warn if some EM's produced NaN.
  if (produced_nan.size() > 0) {
    std::cout << std::endl;
    std::cout << "WARNING: Some iterations produced NaN." << std::endl;
    std::list<long>::iterator it;
    for (it = produced_nan.begin(); it != produced_nan.end(); it++) {
      std::cout << *it << ", ";
    }
    std::cout << std::endl;
  }
}

// below we test the fit the fit of the theoretical distribution to the data, i.e. we simulate data with all visible and check it per node,
// that is what was used to produce the histograms in the paper

void parameter_cloud(Tree &T, Model &Mod, long Nrep, long length, double eps, Parameters &Parsim){

  long iter;

  double likel;
   
	
  Parameters Par;
  Alignment align;
  Counts data;

  double eps_pseudo = 0.001;     // Amount added to compute the pseudo-counts.

  // Initialize the parameters for simulation of K81 data for testing
  create_parameters(Par, T);

  // Obtaining the distribution of estimated parameters with EM
  
  std::ofstream estpar;
  estpar.open("est-par.dat", std::ios::out);
  estpar.precision(15);

  std::vector<double> param;
  for (iter=0; iter < Nrep; iter++) {
    random_data(T, Mod, Parsim, length, align);
    get_counts(align, data);
    add_pseudocounts(eps_pseudo, data);

    // fixed initial conditions
    initial_parameters(Mod, Par);

    // Runs EM
    likel = EMalgorithm(T, Mod, Par, data, eps);

    // Choses the best permutation.
    guess_permutation(T, Mod, Par);

    get_free_param_vector(T, Mod, Par, param);

    for (unsigned long k=0; k < param.size(); k++) {
      estpar << param[k] << "  ";
    }
    estpar << std::endl;
  }

}




