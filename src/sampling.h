/*
 *  sampling.h
 *  
 *  Created by Ania M. Kedzierska on 11/11/11.
 *  Copyright 2011 Politecnic University of Catalonia, Center for Genomic Regulation.  This is program can be redistributed, modified or else as given by the terms of the GNU General Public License. 
 *  
 */

#ifndef __SAMPLING_H__
#define __SAMPLING_H__

#include "tree.h"
#include "alignment.h"
#include "parameters.h"
#include "model.h"
#include "matrix.h"

// long multinomial(Root &p, long nalpha);

void random_parameters(Model &Mod, Parameters &Par);
void random_parameters_length(Tree &T, Model &Mod, Parameters &Par);

void random_fake_counts(Tree &T, double N, Counts &data, Parameters &Parsim);
void random_data(Tree &T, Model &Mod, Parameters &Par, long length, Alignment &align);
void random_full_counts(Tree &T, Model &Mod, Parameters &Par, long length, Matrix &F);

bool check_model_matrix(Model &Mod, TMatrix &tm);
bool check_DLC_matrix(TMatrix &tm);
bool check_matrix_length(double len, TMatrix &tm);
bool check_matrix_stochastic(TMatrix &tm);

#endif
