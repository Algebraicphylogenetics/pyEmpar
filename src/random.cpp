/*
 *  random.cpp
 *  
 *  Created by Ania M. Kedzierska on 11/11/11.
 *  Copyright 2011 Politecnic University of Catalonia, Center for Genomic Regulation.  This is program can be redistributed, modified or else as given by the terms of the GNU General Public License. 
 *  
 */

#include <ctime>

#include "random.h"

#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_real_distribution.hpp>
#include <boost/random/uniform_int_distribution.hpp>
#include <boost/random/discrete_distribution.hpp>


boost::random::mt19937 random_gen((unsigned int) std::time(0));  // Random number generator

boost::random::uniform_real_distribution<double> uni_real_dist;  // Uniform distribution
boost::random::uniform_int_distribution<long> uni_int_dist;      // Uniform distribution
boost::random::discrete_distribution<long, double> disc_dist;    // Discrete distribution


// Initialize the seed of the random number generator
void random_initialize(unsigned int seed) {
  random_gen.seed(seed);
}

// Initialize with a time derived seed
void random_initialize() {
  random_initialize((unsigned int) std::time(0));
}


// uniformly random real in [a,b)
double uniform_real(double a, double b) {
  boost::random::uniform_real_distribution<double>::param_type param(a, b);
  uni_real_dist.param(param);
  return uni_real_dist(random_gen);
}

// random integer in[a,b]
  long uniform_int(long a, long b) {
  boost::random::uniform_int_distribution<long>::param_type param(a, b);
  uni_int_dist.param(param);
  return uni_int_dist(random_gen);
}

// random integer between 0 and p.size()-1 according to the probabilities in p.
long discrete(std::vector<double> &p) {
  boost::random::discrete_distribution<long, double>::param_type param(p);
  disc_dist.param(param);
  return disc_dist(random_gen);
}
