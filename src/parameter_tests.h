/*
 *  parameter_test.h
 *  
 *  Created by Ania M. Kedzierska on 11/11/11.
 *  Copyright 2011 Politecnic University of Catalonia, Center for Genomic Regulation.  This is program can be redistributed, modified or else as given by the terms of the GNU General Public License. 
 *  
 */
#ifndef __PARAMETER_TESTS_H__
#define __PARAMETER_TESTS_H__

#include "tree.h"
#include "model.h"
#include <vector>
#include <string>


enum ParameterTestType {PT_KL, PT_MCEXACT};

void parameter_test(Tree &T, Model &Mod, long Nrep, long length, double eps, std::vector<double> &pvals, std::string data_prefix , bool save_mc_exact);

void parameter_cloud(Tree &T, Model &Mod, long Nrep, long length, double eps, Parameters &Parsim);


struct ci_binom {
	double lower;
	double upper;
};

#endif
