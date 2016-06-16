/*
 *  permutation.h
 *  
 *  Created by Ania M. Kedzierska on 11/11/11.
 *  Copyright 2011 Politecnic University of Catalonia, Center for Genomic Regulation.  This is program can be redistributed, modified or else as given by the terms of the GNU General Public License. 
 *  
 */


#ifndef __PERMUTATION_H__
#define __PERMUTATION_H__

#include <list>
#include <vector>

#include "parameters.h"
#include "tree.h"

typedef std::vector<long> Permutation;
struct Model;

void permute_at_node(Tree &T, long n, Parameters &Par, Permutation &perm);
void guess_permutation(Tree &T, Model &Mod, Parameters &Par);

void permute_rows(TMatrix &tm, Permutation &perm);
void permute_cols(TMatrix &tm, Permutation &perm);
void permute_root(Root &r, Permutation &perm);

bool is_permutation(Permutation &perm);
long permutation_sign(Permutation &perm);

#endif
