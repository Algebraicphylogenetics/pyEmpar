/*
 *  model.h
 *  
 *  Created by Ania M. Kedzierska on 11/11/11.
 *  Copyright 2011 Politecnic University of Catalonia, Center for Genomic Regulation.  This is program can be redistributed, modified or else as given by the terms of the GNU General Public License. 
 *  
 */


#ifndef __MODEL_H__
#define __MODEL_H__

#include <string>
#include <list>

#include "model_k81.h"
#include "model_k80.h"
#include "model_ssm.h"
#include "model_gmm.h"
#include "model_jc.h"

#include "parameters.h"
#include "permutation.h"


enum ModelType {JC, K80, K81, SSM, GMM};

// Struct that stores pointers to the model specific things.
// This struct should be the interface for the model specific functions you use !
// If you want to add a model xxx, create model_xxx.h and model_xxx.cpp, add #include "model_xxx.h"
// here, and modify the function create_model to take the new model into account.

struct Model {
  long nalpha;                                         // Number of bases
  std::string name;                                    // Name of the model
  ModelType m;                                         // ID of the model
  long rdf;                                            // Root degrees of freedom
  long df;                                             // Degrees of freedom
  long np;                                             // Number of parameters
  void (*mle_edge)(TMatrix&, TMatrix&);                // MLE for an edge
  void (*mle_root)(Root&, Root&);                      // MLE for root
  void (*random_edge)(TMatrix&);                       // Random edge
  void (*random_edge_length)(double, TMatrix&);        // Random edge with fixed length
  void (*random_edge_bio_length)(double, TMatrix&);    // Random biol. meaningful edge with fixed length
  void (*random_root)(Root&);                          // Random root
  long (*matrix_structure)(long, long);                // Describes the matrix structure
  long (*root_structure)(long);                        // Describes the root distr. structure.
  void (*list_permutations)(std::list<Permutation>&);  // Lists all compatible permutations.
};


Model create_model(std::string model);



#endif
