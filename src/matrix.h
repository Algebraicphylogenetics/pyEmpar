/*
 *  matrix.h
 *  
 *  Created by Ania M. Kedzierska on 11/11/11.
 *  Copyright 2011 Politecnic University of Catalonia, Center for Genomic Regulation.  This is program can be redistributed, modified or else as given by the terms of the GNU General Public License. 
 *  
 */

#ifndef __MATRIX_H__
#define __MATRIX_H__

#include <fstream>
#include <iostream>
#include <vector>

// This struct encodes the data for a matrix. It needs to be created and destroyed manually.
struct Matrix {
  long nrows;            // Number of rows.
  long ncols;            // number of columns.
  double **m;            // Data in a doubly indexed array: m[row][col].
};



Matrix create_matrix(long nrows, long ncols);
void delete_matrix(Matrix &mat);

void print_matrix(Matrix &M, long rows, long cols, std::ostream &s = std::cout);
void save_matrix(Matrix &M, long rows, long cols, std::string const &fname);

void print_matrix(std::vector<std::vector<double> > &M, std::ostream &s);
void save_matrix(std::vector<std::vector<double> > &M, std::string const &fname);

#endif
