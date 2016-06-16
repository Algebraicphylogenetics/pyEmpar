/*
 *  state.h
 *  
 *  Created by Ania M. Kedzierska on 11/11/11.
 *  Copyright 2011 Politecnic University of Catalonia, Center for Genomic Regulation.  This is program can be redistributed, modified or else as given by the terms of the GNU General Public License. 
 *  
 */

#ifndef __STATE_H__
#define __STATE_H__

#include <vector>
#include <string>

// This struct encodes a state on a set of nodes (len is the number of nodes).
// The state at every node is given by an integer between 0 and nalpha-1.

// To every state, we assign an index (integer). It is useful in loops, e.g. len=3 and nalpha =4, the assignment goes as follows:
//
// 0 ---> (0,0,0)
// 1 ---> (0,0,1)
// 4 ---> (0,1,0)
// 7 ---> (0,1,3)
//     ...
//
// the index associated to a state is an integer between 0 and nalpha^len - 1.

struct State {
  std::vector<long> s;   // vector of states
  long len;              // Length of the array (i.e. number of nodes)
  long nalpha;           // number of states per position.
};


void create_state(State &sta, long length, long nalpha);
void string2state(std::string &str, State &sta);
void state2string(State &sta, std::string &str);
void index2state(long id, State &sta);
long state2index(State &sta);
void permute_state(std::vector<long> &d, State &sta1, State &sta2);

#endif
