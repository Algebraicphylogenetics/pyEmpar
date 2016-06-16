/*
 *  parameters.h
 *  
 *  Created by Ania M. Kedzierska on 11/11/11.
 *  Copyright 2011 Politecnic University of Catalonia, Center for Genomic Regulation.  This is program can be redistributed, modified or else as given by the terms of the GNU General Public License. 
 *  
 */


#ifndef __PARAMETERS_H__
#define __PARAMETERS_H__

#include <string>
#include <iostream>
#include <ostream>

#include "tree.h"

typedef std::vector<std::vector<double> > TMatrix;
typedef std::vector<double> Root;


// This struct stores the parameters for a given tree topology. 
//It automatically gets destroyed when not needed
struct Parameters {
  long nedges;              // Number of edges.
  long nalpha;              // Number of letters in the alphabet.
  std::vector<TMatrix> tm; // Array of transition matrices. Triply indexed: tm[edge][row][col].
  Root r;                   // Root distribution. An array of length nalpha.
};




void create_parameters(Parameters &Par, Tree &T);

double branch_length(TMatrix &tm, long nalpha);
void branch_lengths(Parameters &Par, std::vector<double> &br);
double parameters_distance(Parameters &Par1, Parameters &Par2);
double parameters_distance_edge(Parameters &Par1, Parameters &Par2, long e);
double parameters_distance_root(Parameters &Par1, Parameters &Par2);

void print_parameters(Parameters &Par, std::ostream &s=std::cout);
void save_parameters(Parameters &Par, std::string const &fname);
void read_parameters(Parameters &Par, std::string const &fname);
double branch_length_error_bound(TMatrix &tm_oryg, TMatrix &tm_est);
double branch_length_error_bound_mult(TMatrix &tm_oryg, TMatrix &tm_est);
double log_multinomial_evaluate(std::vector<double> &p, std::vector<double> &data, long N);
double get_mult(std::vector<double> &tm_oryg, std::vector<double> &counts_est, long N);
double bound_mult(TMatrix &tm_oryg, double pr_mult, long N);

void copy_parameters(Parameters &source, Parameters &target);

#endif
