/*
 *  fisher.h
 *  
 *  Created by Ania M. Kedzierska on 11/11/11.
 *  Copyright 2011 Politecnic University of Catalonia, Center for Genomic Regulation.  This is program can be redistributed, modified or else as given by the terms of the GNU General Public License. 
 *  
 */


#ifndef __FISHER_H__
#define __FISHER_H__

#include "tree.h"
#include "model.h"
#include "parameters.h"
#include "matrix.h"
#include "alignment.h"

typedef std::vector<double> Array1;                              // 1-dimensional array
typedef std::vector<std::vector<double> > Array2;                // 2-dimensional array
typedef std::vector<std::vector<std::vector<double> > > Array3;  // 3-dimensional array

void Fisher_information(Tree &T, Model &Mod, Parameters &Par, long N, Matrix &I);
void Observed_Fisher_information(Tree &T, Model &Mod, Parameters &Par, Counts &data, Array2 &I);

// void Fisher_information_approx(Tree &T, Model &Mod, Parameters &Par, long N, Matrix &I);

void full_MLE_covariance_matrix(Tree &T, Model &Mod, Parameters &Par, long N, Array2 &Cov);
void full_MLE_observed_covariance_matrix(Tree &T, Model &Mod, Parameters &Par, Counts &data, Array2 &Cov);

void branch_covariance_matrix(Model &Mod, Array2 &Covfull, long b, Array2 &Covbr);
void branch_inverted_covariance_matrix(Model &Mod, Array2 &Covfull, long b, Array2 &Covbri);

void get_param_vector(Tree &T, Model &Mod, Parameters &Par, std::vector<double> &param);
void get_free_param_vector(Tree &T, Model &Mod, Parameters &Par, std::vector<double> &param);
void get_branch_free_param_vector(Tree &T, Model &Mod, Parameters &Par, long br, std::vector<double> &param);


#endif
