/*
 *  model.cpp
 *  
 *  Created by Ania M. Kedzierska on 11/11/11.
 *  Copyright 2011 Politecnic University of Catalonia, Center for Genomic Regulation.  This is program can be redistributed, modified or else as given by the terms of the GNU General Public License. 
 *  
 */


#include "model.h"



// Takes a model name as input, and fills the struct Model with pointers
// to the right model specific functions.
// If we want to add a new model, this function must be modified.

Model create_model(std::string model) {
  Model M;
  if (model == "k81") {
    M.name = "Kimura 81";
    M.m = K81;
    M.df = 3;
    M.np = 4;
    M.rdf = 0;
    M.mle_root = &K81_mle_root;
    M.mle_edge = &K81_mle_edge;
    M.random_root = &K81_random_root;
    M.random_edge = &K81_random_edge;
    M.random_edge_length = &K81_random_edge_length;
    M.random_edge_bio_length = &K81_random_edge_bio_length;
    M.matrix_structure = &K81_matrix_structure;
    M.root_structure = &K81_root_structure;
    M.list_permutations = &K81_list_permutations;
    M.nalpha = 4;
  } else if (model == "k80") {
    M.name = "Kimura 80";
    M.m = K80;
    M.df = 2;
    M.np = 3;
    M.rdf = 0;
    M.mle_root = &K80_mle_root;
    M.mle_edge = &K80_mle_edge;
    M.random_root = &K80_random_root;
    M.random_edge = &K80_random_edge;
    M.random_edge_length = &K80_random_edge_length;
    M.random_edge_bio_length = &K80_random_edge_bio_length;
    M.matrix_structure = &K80_matrix_structure;
    M.root_structure = &K80_root_structure;
    M.list_permutations = &K80_list_permutations;
    M.nalpha = 4;
  } else if (model == "ssm") {
    M.name = "Strand Symmetric";
    M.m = SSM;
    M.df = 6;
    M.np = 8;
    M.rdf = 1;
    M.mle_root = &SSM_mle_root;
    M.mle_edge = &SSM_mle_edge;
    M.random_root = &SSM_random_root;
    M.random_edge = &SSM_random_edge;
    M.random_edge_length = &SSM_random_edge_length;
    M.random_edge_bio_length = &SSM_random_edge_bio_length;
    M.matrix_structure = &SSM_matrix_structure;
    M.root_structure = &SSM_root_structure;
    M.list_permutations = &SSM_list_permutations;
    M.nalpha = 4;
  } else if (model == "jc") {
    M.name = "Jukes Cantor";
    M.m = JC;
    M.df = 1;
    M.np = 2;
    M.rdf = 0;
    M.mle_root = &JC_mle_root;
    M.mle_edge = &JC_mle_edge;
    M.random_root = &JC_random_root;
    M.random_edge = &JC_random_edge;
    M.random_edge_length = &JC_random_edge_length;
    M.random_edge_bio_length = &JC_random_edge_bio_length;
    M.matrix_structure = &JC_matrix_structure;
    M.root_structure = &JC_root_structure;
    M.list_permutations = &JC_list_permutations;
    M.nalpha = 4;
  } else if (model == "gmm") {
    M.name = "General Markov";
    M.m = GMM;
    M.df = 12;
    M.np = 16;
    M.rdf = 3;
    M.mle_root = &GMM_mle_root;
    M.mle_edge = &GMM_mle_edge;
    M.random_root = &GMM_random_root;
    M.random_edge = &GMM_random_edge;
    M.random_edge_length = &GMM_random_edge_length;
    M.random_edge_bio_length = &GMM_random_edge_bio_length;
    M.matrix_structure = &GMM_matrix_structure;
    M.root_structure = &GMM_root_structure;
    M.list_permutations = &GMM_list_permutations;
    M.nalpha = 4;
  }
  return M;
}

