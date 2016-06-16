/*
 *  tree.cpp
 *  
 *  Created by Ania M. Kedzierska on 11/11/11.
 *  Copyright 2011 Politecnic University of Catalonia, Center for Genomic Regulation.  This is program can be redistributed, modified or else as given by the terms of the GNU General Public License. 
 *  
 */


#include <iostream>
#include <fstream>
#include <list>
#include <algorithm>

#include "Newickform.h"
#include "seqUtil.h"
#include "tree.h"




// Stores a pointer to every node in an array. First the leaves, then the rest.
void flatten_tree(newick_node* root, std::list<newick_node*> &Lleafs, std::list<newick_node*> &Lhidden) {
  newick_child* ch;

  if (root == NULL) {
    return;
  }

  ch = root->child;
  if (ch == NULL) {
    Lleafs.push_back(root);
  } else {
    Lhidden.push_back(root);
  }

  while (ch != NULL) {
    flatten_tree(ch->node, Lleafs, Lhidden);
    ch = ch->next;
  }  
}


// Fills a list of all edges in the tree, from lists of newick_node pointers.
void list_edges(std::vector<newick_node*> &nodes, std::list<Edge> &edges) {
  unsigned int i,j;
  newick_child* ch;
  Edge e;

  for (i=0; i < nodes.size(); i++) {
    ch = nodes[i]->child;
    while (ch != NULL) {
      for(j=0; j < nodes.size(); j++) {
        if (nodes[j] == ch->node) break;
      }
      e.s = i;
      e.t = j;
	  e.br = 0;
      edges.push_back(e);
      ch = ch->next;
    }
  }  
}


// Reads a Newick tree from a file and parses it, returning a Tree struct.
void read_tree(Tree &T, std::string fname, long nalpha) {
  long i;
   newick_node *root;
  //  FILE *f;

  // Memory initialization for the Newick code.
  seqMemInit();

  std::ifstream ftree;
  std::string tree_string, line;

  ftree.open(fname.c_str(), std::fstream::in);
  tree_string = "";
  if (!ftree.is_open()) {
    std::cout << "Cannot open the tree file." << std::endl;
    exit(1);
  } 
  while (!ftree.eof()) {
    // read a line
    ftree >> line;
    // trim spaces at the begining and the end
    line.erase(remove_if(line.begin(), line.end(), ::isspace), line.end());
    tree_string.append(line);
  }
  ftree.close();

  // copies the C++ string into a C-style string for the parser 
  char *cstr = new char[tree_string.size()+1];
   std::copy(tree_string.begin(), tree_string.end(), cstr);
   cstr[tree_string.size()] = '\0';

  // parses the string
  root = parseTree(cstr); 

  // releases memory for the C string.
  delete[] cstr;

  // Putting pointers to the nodes in a vector. First leaves, then hidden nodes.
  std::list<newick_node*> list_leafs;
  std::list<newick_node*> list_hidden;
  std::vector<newick_node*> nodes;
  std::list<Edge> edges;

  flatten_tree(root, list_leafs, list_hidden);
  nodes.insert(nodes.begin(), list_leafs.begin(), list_leafs.end());
  nodes.insert(nodes.end(), list_hidden.begin(), list_hidden.end());

  // Listing all the edges. Every edge is a struct with source and target.
  // Source and target are the indices of the node in the vector nodes.
  list_edges(nodes, edges);

  // Filling the Tree struct.
  T.edges.resize(0);
  T.edges.insert(T.edges.begin(), edges.begin(), edges.end());
  T.nalpha = nalpha;
  T.nedges = T.edges.size();
  T.nleaves = list_leafs.size();
  T.nhidden = list_hidden.size();
  T.nnodes = T.nleaves + T.nhidden;

  T.nsthidden = 1;
  for(i=0; i < T.nhidden; i++) {
    T.nsthidden = T.nsthidden*T.nalpha;
  }

  T.nstleaves = 1;
  for(i=0; i < T.nleaves; i++) {
    T.nstleaves = T.nstleaves*T.nalpha;
  }


  for (i=0; i < (long) T.edges.size(); i++) {
    T.edges[i].br = (double) ((nodes[T.edges[i].t])->dist);
  }

  T.names.resize(nodes.size());
  for (i=0; i < T.nleaves; i++) {
    ((T.names)[i]).assign((nodes[i])->taxon);
  }

  T.tree_name = fname;

  create_state(T.shidden, T.nhidden, T.nalpha);
  create_state(T.sleaves, T.nleaves, T.nalpha);

  // Freeing memory for the Newick code.
  seqFreeAll();
}


// Computes the valence of a node: the number of incoming + outgoing edges.
long valence(Tree &T, long node) {
  long e;
  long val=0;
  for(e=0; e < T.nedges; e++) {
    if (T.edges[e].s == node || T.edges[e].t == node){
      val = val + 1;
    }
  }
  return val;
}


// Prints the edges of a tree.
void print_tree(Tree &T) {
  unsigned int i;
  std::cout << "nodes:   " << T.nnodes << std::endl;
  std::cout << "nleaves: " << T.nleaves << std::endl;
  std::cout << "nedges:  " << T.nedges << std::endl;
  std::cout << "Edges:" << std::endl;
  for (i=0; i < T.edges.size(); i++) {
    std::cout << "  (" << T.edges[i].s << ", " << T.edges[i].t << ")";
    std::cout << "  " << T.edges[i].br << std::endl;
  }
  std::cout << "Node names: ";
  for (i=0; i < T.nleaves; i++) {
    std::cout << T.names[i] << " ";
  }
  std::cout << std::endl;
  std::cout << std::endl;
}


// Lists all outgoing edges from a given node.
void list_outgoing_edges(Tree &T, long node, std::list<long> &L) {
  L.clear();
  for(long e=0; e < T.nedges; e++) {
    if (T.edges[e].s == node) {
      L.push_back(e);  
    }
  }
}


// recursively prints the newick tree from the root.
void print_newick_tree_rec(Tree &T, std::vector<double> br, long root) {
  long e;
  bool first;
  if (root < T.nleaves) {   // if root is already a leaf
    std::cout << T.names[root];
  } else {                  // otherwise print all the subtrees
    std::cout << "(";
    first = true;
    for (e=0; e < T.nedges; e++) {
      if (T.edges[e].s == root) {
        if (!first) {
          std::cout << ",";    // If not the first child, print a comma.
        }
        print_newick_tree_rec(T, br, T.edges[e].t);
	std::cout << ":" << br[e];
        first = false;
      } 
    }
    std::cout << ")";
  }
}



// Prints the tree T in Newick format, replacing the branch lengths by
// the ones in br.
void print_newick_tree(Tree &T, std::vector<double> br) {
  print_newick_tree_rec(T, br, T.nleaves);  // The root is the node nleaves.
  std::cout << std::endl << std::endl;
}
