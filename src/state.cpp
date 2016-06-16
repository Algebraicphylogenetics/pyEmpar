/*
 *  state.cpp
 *  
 *  Created by Ania M. Kedzierska on 11/11/11.
 *  Copyright 2011 Politecnic University of Catalonia, Center for Genomic Regulation.  This is program can be redistributed, modified or else as given by the terms of the GNU General Public License. 
 *  
 */

#include <iostream>
#include <cstdlib>

#include "state.h"

// creates a state of given length and ginven number of letters.
void create_state(State &sta, long length, long nalpha) {
  sta.len = length;
  sta.nalpha = nalpha;
  sta.s.resize(length);
}

// puts the state associated to an index in sta.
void index2state(long id, State &sta) {
  int i;
  for (i=0; i<sta.len; i++) {
    sta.s[sta.len - i - 1] = id % sta.nalpha; // rest of dividing by nalpha.
    id = (id - sta.s[sta.len - i - 1]) / sta.nalpha;
  }
}

// returns the index associated to a state.
long state2index(State &sta) {
  int i;
  long id = 0;
  for (i=0; i < sta.len; i++) {
    id = sta.nalpha*id + sta.s[i];
  }
  return id;
}

// converts a state into a string.
void state2string(State &sta, std::string &str) {
  char bases[4] = {'a', 'c', 'g', 't'};
  int i;
  if (sta.nalpha == 4) {
    str.resize(0);
    for (i=0; i < sta.len; i++) {
      str.push_back(bases[sta.s[i]]);
    }
  } else {
    std::cout << "Can't convert a state with more than 4 letters to a string." << std::endl;
    exit(1);
  }
}

// converts a string into a state.
void string2state(std::string &str, State &sta) {
  unsigned int i;
  char c;
  sta.nalpha = 4;
  sta.len = str.size();
  sta.s.resize(sta.len);
  for (i=0; i < str.size(); i++) {
    c = tolower(str[i]);
    if (c == 'a') sta.s[i] = 0;
    if (c == 'c') sta.s[i] = 1;
    if (c == 'g') sta.s[i] = 2;
    if (c == 't') sta.s[i] = 3;
  }
}

// Permutes sta1 and stores into sta2, such that sta2.s[d[i]] = sta1.s[i]
void permute_state(std::vector<long> &d, State &sta1, State &sta2) {
  long i;
  if ((long) d.size() != sta1.len || sta1.len != sta2.len) {
    std::cout << "Permutation and states are not compatible !" << std::endl;
    exit(-1);
  }
  for(i=0; i < sta1.len; i++) {
    sta2.s[d[i]] = sta1.s[i];
  } 
}

