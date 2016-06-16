/*
 * model_k81.h
 *  
 *  Created by Ania M. Kedzierska on 11/11/11.
 *  Copyright 2011 Politecnic University of Catalonia, Center for Genomic Regulation.  This is program can be redistributed, modified or else as given by the terms of the GNU General Public License. 
 *  
 */

#ifndef __MODEL_K81_H__
#define __MODEL_K81_H__

#include "parameters.h"
#include "permutation.h"
#include <list>

void K81_mle_edge(TMatrix &N, TMatrix &tm);
void K81_mle_root(Root &s, Root &r);
void K81_random_edge(TMatrix &tm);
void K81_random_edge_length(double len, TMatrix &tm);
void K81_random_edge_bio_length(double len, TMatrix &tm);
void K81_random_root(Root &r);
long K81_matrix_structure(long i, long j);
long K81_root_structure(long i);
void K81_list_permutations(std::list<Permutation> &L);

#endif
