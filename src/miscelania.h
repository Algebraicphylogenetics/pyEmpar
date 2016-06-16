/*
 *  miscelania.h
 *  
 *  Created by Ania M. Kedzierska on 11/11/11.
 *  Copyright 2011 Politecnic University of Catalonia, Center for Genomic Regulation.  This is program can be redistributed, modified or else as given by the terms of the GNU General Public License. 
 *  
 */

#ifndef __MISCELANIA_H__
#define __MISCELANIA_H__

#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include "parameters.h"

std::string strip_extension(std::string fname);
std::string make_parameters_filename(std::string fastafile);

void print_vector(std::vector<double> &v, std::ostream &s = std::cout);
void save_vector(std::vector<double> &v, std::string const &fname);

void swap(double &a, double &b);
long max_in_col(TMatrix &tm, long col);

#endif
