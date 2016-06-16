/*
 *  em.h
 *  
 *  Created by Ania M. Kedzierska on 11/11/11.
 *  Copyright 2011 Politecnic University of Catalonia, Center for Genomic Regulation.  This is program can be redistributed, modified or else as given by the terms of the GNU General Public License. 
 *  
 */

#ifndef __EM_H__
#define __EM_H__

#include <vector>

#include "tree.h"
#include "parameters.h"
#include "alignment.h"
#include "model.h"
#include "state_list.h"
#include "matrix.h"


double joint_prob(Tree &T, Parameters &Par, long idleaf, long idhid);
double joint_prob_leaves(Tree &T, Parameters &Par, long idleaf);
double log_likelihood(Tree &T, Parameters &Par, Counts &data);
double KL_divergence(Tree &T, Parameters &Par1, Parameters &Par2);
double entropy(Tree &T, Parameters &Par);

double joint_prob_fast(Tree &T, Parameters &Par, StateList &sl, long idleaf, long idhid);
double joint_prob_leaves_fast(Tree &T, Parameters &Par, StateList &sl, long idleaf);
double log_likelihood_fast(Tree &T, Parameters &Par, Counts &data, StateList &sl);
double KL_divergence_fast(Tree &T, Parameters &Par1, Parameters &Par2, StateList &sl);
double entropy_fast(Tree &T, Parameters &Par, StateList &sl);

void two_node_marginalization(Tree &T, Matrix &F, long a, long b, StateList &sl, TMatrix &N);
void one_node_marginalization(Tree &T, Matrix &F, long a, StateList &sl, Root &s);

double EMalgorithm(Tree &T, Model &Mod, Parameters &Par, Counts &data, double eps, bool silent=false);
void MLE_all_obs(Tree &T, Model &Mod, Parameters &Par, Matrix &F, StateList &sl);

void initial_parameters(Model &Mod, Parameters &Par);

#endif
