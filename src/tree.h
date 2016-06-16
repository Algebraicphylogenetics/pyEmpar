/*
 *  tree.h
 *  
 *  Created by Ania M. Kedzierska on 11/11/11.
 *  Copyright 2011 Politecnic University of Catalonia, Center for Genomic Regulation.  This is program can be redistributed, modified or else as given by the terms of the GNU General Public License. 
 *  
 */

#ifndef __TREE_H__
#define __TREE_H__

#include <string>
#include <vector>
#include <list>
#include "state.h"

// Stores an edge. s is the source node, t the target and br a branch length.
struct Edge {
  long s;                            // source node
  long t;                            // target node
  double br;                         // branch length
};


// Struct that stores a tree as a list of edges.
// The auxiliar states shidden and sleaves are used internally on the EM algorithm.
struct Tree { 
  long nalpha;                       // number of states per node
  long nleaves;                      // number of leaves
  long nnodes;                       // number of nodes
  long nhidden;                      // number of hidden nodes (nnodes - nleaves)
  long nstleaves;                    // number of states on the leaves: nalpha^nleaves
  long nsthidden;                    // number of hidden states: nalpha^nhidden
  long nedges;                       // number of edges
  std::vector<Edge> edges;           // vector with the edges.
  std::vector<std::string> names;    // node names
  std::string tree_name;             // filename from where the tree was read.
  
  // Auxiliar data for using inside the em algorithm.
  State shidden;                     // stores a state on the hidden nodes
  State sleaves;                     // Stores a state on the leaves
};


void read_tree(Tree &T, std::string fname, long nalpha=4);
void print_tree(Tree &T);
void print_newick_tree(Tree &T, std::vector<double> br);

void list_outgoing_edges(Tree &T, long node, std::list<long> &L);
long valence(Tree &T, long node);


#endif
