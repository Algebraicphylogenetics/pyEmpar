/*
 * state_list.h
 *  
 *  Created by Ania M. Kedzierska on 11/11/11.
 *  Copyright 2011 Politecnic University of Catalonia, Center for Genomic Regulation.  This is program can be redistributed, modified or else as given by the terms of the GNU General Public License. 
 *  
 */

#ifndef __STATE_LIST_H__
#define __STATE_LIST_H__


#include <vector>
#include "state.h"
#include "tree.h"


// The struct StateList encodes a list of all states on the hidden nodes and leaves
// of a tree. It is used internally in the EM algorithm for speeding up the conversion
// from a state index to a vector of states on the nodes.
// e.g. h[i][k] gives the state in which the k-th hidden node is, when the global
// state on hidden nodes is i.

struct StateList {
  long nalpha;                       // number of states per node
  long nhidden;                      // number of hidden nodes
  long nleaves;                      // number of leaf nodes
  long nsthidden;                    // number of states on hidden nodes
  long nstleaves;                    // number of states on leaves
  std::vector<std::vector<int> > h;  // array of states on hidden nodes.
  std::vector<std::vector<int> > l;  // array of states on leaves.
};


void create_state_list(StateList &sl, Tree &T);

#endif
