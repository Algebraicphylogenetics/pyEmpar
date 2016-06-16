/*
 *  miscalenia.cpp
 *  
 *  Created by Ania M. Kedzierska on 11/11/11.
 *  Copyright 2011 Politecnic University of Catalonia, Center for Genomic Regulation.  This is program can be redistributed, modified or else as given by the terms of the GNU General Public License. 
 *  
 */

#include <cstdlib>
#include <iostream>

#include "miscelania.h"

// cut-off the extension from the input filename.
std::string strip_extension(std::string fname) {
  int i;
  std::string ret;
  for(i=fname.size()-1; i>= 0; i--) {
    if (fname[i] == '.') break;
  }
  if (i<= 0) ret = fname;
  else ret.assign(fname, 0, i);
  return ret;
}

// Create the ame of the output file, replacing ".fa" by ".dat".
std::string make_parameters_filename(std::string fastafile) {
  return strip_extension(fastafile) + ".dat";
}


// prints a vector of doubles in a row. Used for printing branch lengths
void print_vector(std::vector<double> &v, std::ostream &s) {
  unsigned int i;
  for(i=0; i < v.size(); i++) {
    s << v[i] << "  ";
  }
  s << std::endl;
  s << std::endl;
}


void save_vector(std::vector<double> &v, std::string const &fname) {
  std::fstream st;
  st.precision(15);

  st.open(fname.c_str(), std::ios::out);
  if (!st.is_open()) {
    std::cout << "Could not open file: " << fname << std::endl;
  }
  print_vector(v, st);
  st.close();
}


void swap(double &a, double &b) {
  double c;
  c = b;
  b = a;
  a = c;
}


long max_in_col(TMatrix &tm, long col) {
  long i0;
  double max;

  i0=col;
  max = tm[col][col];
  for (unsigned long i=0; i < tm.size(); i++) {
    if (tm[i][col] > max) {
      i0 = i;
      max = tm[i][col];
    }
  }
  return i0;
}
