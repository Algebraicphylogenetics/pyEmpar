/*
 *  matrix.cpp
 *  
 *  Created by Ania M. Kedzierska on 11/11/11.
 *  Copyright 2011 Politecnic University of Catalonia, Center for Genomic Regulation.  This is program can be redistributed, modified or else as given by the terms of the GNU General Public License. 
 *  
 */

#include <iostream>
#include "matrix.h"


// Assigns memory for a matrix of given number of rows and columns
Matrix create_matrix(long nrows, long ncols) {
  Matrix mat;
  int i;

  mat.nrows = nrows;
  mat.ncols = ncols;
  mat.m = new double*[nrows];

  for(i=0; i<nrows; i++) {
    mat.m[i] = new double[nrows];
  }
  return mat;
}

// Frees the memory assigned to the matrix
void delete_matrix(Matrix &mat) {
  int i;

  for(i=0; i < mat.nrows; i++) {
    delete[] mat.m[i];
  }
  delete[] mat.m;
}




// Prints a matrix
void print_matrix(Matrix &M, long rows, long cols, std::ostream &s) {
  long i,j;
  for(i=0; i < rows; i++) {
    for(j=0; j < cols; j++) {
      s << M.m[i][j] << " ";
    }
    s << std::endl;
  }
}

// Saves a matrix into a file
void save_matrix(Matrix &M, long rows, long cols, std::string const &fname) {
  std::fstream st;
  st.precision(15);

  st.open(fname.c_str(), std::ios::out);
  if (!st.is_open()) {
    std::cout << "Could not open file: " << fname << std::endl;
  }
  print_matrix(M, rows, cols, st);
  st.close();
}



// Prints a matrix
void print_matrix(std::vector<std::vector<double> > &M, std::ostream &s) {
  unsigned long i,j;
  for(i=0; i < M.size(); i++) {
    for(j=0; j < M[i].size(); j++) {
      s << M[i][j] << " ";
    }
    s << std::endl;
  }
}

// Saves a matrix into a file
void save_matrix(std::vector<std::vector<double> > &M, std::string const &fname) {
  std::fstream st;
  st.precision(15);

  st.open(fname.c_str(), std::ios::out);
  if (!st.is_open()) {
    std::cout << "Could not open file: " << fname << std::endl;
  }
  print_matrix(M, st);
  st.close();
}

