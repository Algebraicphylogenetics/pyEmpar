/*
 *  parameters.cpp
 *  
 *  Created by Ania M. Kedzierska on 11/11/11.
 *  Copyright 2011 Politecnic University of Catalonia, Center for Genomic Regulation.  This is program can be redistributed, modified or else as given by the terms of the GNU General Public License. 
 *  
 */


#include "state_list.h"



// Assigns memory and fills the struct StateList from the information in T.
// After this, sl.l containes the states on leaves and sl.h the states on hidden
// nodes. Both states are obtained from the index via the index2state function in state.h.
void create_state_list(StateList &sl, Tree &T) {
  long i, j;
  State sta;

  sl.nalpha = T.nalpha;


  // Initialize the number of hidden states
  sl.nhidden = T.nhidden;
  sl.nsthidden = 1;
  for (i=0; i < sl.nhidden; i++) {
    sl.nsthidden = sl.nsthidden * sl.nalpha;
  }

  // fill in the state list for the hidden nodes.
  sl.h.resize(sl.nsthidden);
  create_state(sta, sl.nhidden, sl.nalpha);
  for (i=0; i < sl.nsthidden; i++) {
    index2state(i, sta);
    sl.h[i].resize(sl.nhidden);
    for(j=0; j < sl.nhidden; j++) {
      sl.h[i][j] = sta.s[j];
    }
  }

  // Initialize the number of leaf states
  sl.nleaves = T.nleaves;
  sl.nstleaves = 1;
  for (i=0; i < sl.nleaves; i++) {
    sl.nstleaves = sl.nstleaves * sl.nalpha;
  }

  // fill in the state list for the leaf nodes.
  sl.l.resize(sl.nstleaves);
  create_state(sta, sl.nleaves, sl.nalpha);
  for (i=0; i < sl.nstleaves; i++) {
    index2state(i, sta);
    sl.l[i].resize(sl.nleaves);
    for(j=0; j < sl.nleaves; j++) {
      sl.l[i][j] = sta.s[j];
    }
  }
}
