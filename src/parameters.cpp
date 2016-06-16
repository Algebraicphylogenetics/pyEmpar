/*
 *  parameters.cpp
 *  
 *  Created by Ania M. Kedzierska on 11/11/11.
 *  Copyright 2011 Politecnic University of Catalonia, Center for Genomic Regulation.  This is program can be redistributed, modified or else as given by the terms of the GNU General Public License. 
 *  
 */

// This file contains a veriaty of functions that are not used in the final tests. I tried the exact test, approaximations, KL divergence, 
// different bounds on the number of parameters etc. The code is left inside for to be possibly of use or played with.

#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cmath>

#include "parameters.h"
#include "em.h"
#include "matrix.h"

// Assigns memory for storing the parameters matching the given tree.
void create_parameters(Parameters &Par, Tree &T) {
  int i,j;

  Par.nalpha = T.nalpha;
  Par.nedges = T.nedges;

  Par.r.resize(T.nalpha);
  Par.tm.resize(T.nedges);
  for(i=0; i<T.nedges; i++) {
    (Par.tm[i]).resize(T.nalpha);
    for(j=0; j<T.nalpha; j++) {
      (Par.tm[i][j]).resize(T.nalpha);
    }
  }
}



// Computes the determinant of an nxn matrix stored in a.
double determinant(TMatrix &a, int n)
{
  int i,j,j1,j2;
  double det = 0;
  TMatrix m;

  if (n < 1) { /* Error */

  } else if (n == 1) { /* Shouldn't be used */
    det = a[0][0];
  } else if (n == 2) {
    det = a[0][0] * a[1][1] - a[1][0] * a[0][1];
  } else {
    det = 0;
    for (j1=0;j1<n;j1++) {
      m.resize(n-1);
      for (i=0;i<n-1;i++)
        m[i].resize(n-1);
      for (i=1;i<n;i++) {
        j2 = 0;
        for (j=0;j<n;j++) {
          if (j == j1)
            continue;
          m[i-1][j2] = a[i][j];
          j2++;
        }
      }
      det += pow(-1.0,1.0+j1+1.0) * a[0][j1] * determinant(m,n-1);
    }
  }
  return(det);
}


// compute the cofactor of the element (row,col)
int GetMinor(TMatrix &src, TMatrix &dest, int row, int col, int order)
{
    // indicate which col and row is being copied to dest
    int colCount=0,rowCount=0;

    for(int i = 0; i < order; i++ )
    {
        if( i != row )
        {
            colCount = 0;
            for(int j = 0; j < order; j++ )
            {
                // when j is not the element
                if( j != col )
                {
                    dest[rowCount][colCount] = src[i][j];
                    colCount++;
                }
            }
            rowCount++;
        }
    }
    return 1;
}


// matrix inversion
// the result is put in Y
void matrix_inverse(TMatrix &A, int order, TMatrix &Y){
    // get the determinant of a
    double det = 1.0/determinant(A,order);
   
    TMatrix minor;
    minor.resize(order-1);
    for(int i=0; i < order-1; i++)
        minor[i].resize(order-1);

    for(int j=0;j<order;j++)
    {
        for(int i=0;i<order;i++)
        {
            // get the co-factor (matrix) of A(j,i)
            GetMinor(A,minor,j,i,order);
            Y[i][j] = det*determinant(minor,order-1);
            if( (i+j)%2 == 1)
                Y[i][j] = -Y[i][j];
        }
    }
}

// Computes the matrix product AB.
void matrix_product(TMatrix &A, TMatrix &B, TMatrix &prod) {
  int rows = A.size();
  int cols = B[0].size();

  prod.resize(rows);
  for (int i=0; i < rows; i++) {
    prod[i].resize(cols);
  }  

  for(int i=0; i < rows; i++) {
    for(int j=0; j < cols; j++) {
      prod[i][j] = 0;
      for(unsigned int k=0; k < B.size(); k++) {
        prod[i][j] = prod[i][j] + A[i][k]*B[k][j];
      }
    }
  }
}

// Computes the matrix A-B.
void matrix_substract(TMatrix &A, TMatrix &B, TMatrix &diff) {
	int rows = A.size();
	int cols = B[0].size();
	
	diff.resize(rows);
	for (int i=0; i < rows; i++) {
		diff[i].resize(cols, 0);
	}  
	
	for(int i=0; i < rows; i++) {
		for(int j=0; j < cols; j++) {
		    diff[i][j] = A[i][j] - B[i][j];
		}
	}
}


// computes the induced L1 norm of a matrix, i.e the maximum
// of the L1 norm of the columns.
double induced_L1_norm(TMatrix &A) {
  double norm=0;
  double newnorm;
  int rows= A.size();
  int cols= A[0].size();

  for (int i=0; i < cols; i++) {
    newnorm=0;
    for(int j=0; j < rows; j++) {
      newnorm = newnorm + fabs(A[i][j]);
    }
    if (newnorm > norm) norm = newnorm;
  }
  return norm;
}

// of the L2 norm of the columns.
double induced_L2_norm(TMatrix &A) {
	double norm=0;
	double newnorm;
	int rows= A.size();
	int cols= A[0].size();
	
	for (int i=0; i < cols; i++) {
		newnorm=0;
		for(int j=0; j < rows; j++) {
			newnorm = newnorm + pow(A[i][j], 2);
		}
		if (newnorm > norm) norm = newnorm;
	}
	return sqrt(norm);
}


// Computes an upper bound for the branch length error.
double branch_length_error_bound(TMatrix &tm_oryg, TMatrix &tm_est) {
	int K = tm_oryg.size();
    TMatrix m_inv, m_diff;
	m_inv.resize(K);
	m_diff.resize(K);
	for(int i=0; i < K; i++) {
		m_inv[i].resize(K, 0);
		m_diff[i].resize(K, 0);
	}
	matrix_inverse(tm_oryg, K, m_inv);
	matrix_substract(tm_oryg, tm_est, m_diff);

	
    double cota = 1./4 * K * induced_L1_norm(m_inv) * induced_L1_norm(m_diff);
	std::cout.precision(15);
	return cota;
}


// Computes an upper bound for the branch length error: version 2 for the mult. probabilities
double branch_length_error_bound_mult(TMatrix &tm_oryg, TMatrix &tm_est) {
	int K = tm_oryg.size(); // vector=1st row for ex. - its size/length
	 
    double vec_diffL2 = 0.0; 
	double temp = 0.0;
	TMatrix mp_inv;// inverse of the original transition matrix
	mp_inv.resize(K);
	for(int i=0; i < K; i++) {
		mp_inv[i].resize(K, 0);
		temp = tm_oryg[0][i] - tm_est[0][i]; 
		vec_diffL2 = vec_diffL2 + pow(temp, 2);
	}
	
	matrix_inverse(tm_oryg, K, mp_inv); // inverse of the original matrix

    double cota = induced_L2_norm(mp_inv) * sqrt(vec_diffL2);
	std::cout.precision(15);
	return cota;
}

// Computes an upper bound for the branch length error: version 2 for the mult. probabilities
double bound_mult(TMatrix &tm_oryg, double pr_mult, long N) {
	int K = tm_oryg.size(); // vector=1st row for ex. - its size/length
	TMatrix mp_inv;// inverse of the original transition matrix
	mp_inv.resize(K);
	for(int i=0; i < K; i++) {
		mp_inv[i].resize(K, 0);
		}
	
	matrix_inverse(tm_oryg, K, mp_inv); // inverse of the original matrix
	
	  double cota = induced_L2_norm(mp_inv) * sqrt(-1* pr_mult/((double) N));
	std::cout.precision(15);
	return cota;
}


double log_multinomial_evaluate(std::vector<double> &p, std::vector<double> &data, long N) {
	unsigned long i;
	
	double ret = (double)N*log((double)N); 
	for (i=0; i < p.size(); i++) {
		  if (data[i] > 0) {
			  ret = ret + data[i] * log(p[i]/data[i]);
 	}
	}
    return ret;
}
// multinomial probability

double get_mult(std::vector<double> &tm_oryg, std::vector<double> &counts_est, long N){ // here they are vectors!!
    double x_mult = 0;
	
    x_mult = exp(log_multinomial_evaluate(tm_oryg, counts_est, N)); // vectors! first row
	return x_mult;
}


// Computes the branch length of a matrix
double branch_length(TMatrix &tm, long nalpha) {
  return -0.25*log(determinant(tm, nalpha));
}


// Computes the branch lengths of the matrices in Par and stores them in a vector.
void branch_lengths(Parameters &Par, std::vector<double> &br) {
  int i;
  br.resize(Par.nedges);
  for (i=0; i < Par.nedges; i++) {
    br[i] = -0.25*log(determinant(Par.tm[i], Par.nalpha));
  }
}

// Computes the L2 distance between the parameters of a given edge e.
double parameters_distance_edge(Parameters &Par1, Parameters &Par2, long e) {
  double c, d;
  long j, k;
  
  d = 0;
  for(j=0; j < Par1.nalpha; j++) {
    for(k=0; k < Par1.nalpha; k++) {
      c = Par1.tm[e][j][k] - Par2.tm[e][j][k];
      d = d + c*c;
    }
  }
  return sqrt(d);
}

// Computes the L2 distance between the roots.
double parameters_distance_root(Parameters &Par1, Parameters &Par2) {
  double d, c;
  long i;
  d = 0;
  for (i=0; i < Par1.nalpha; i++) {
    c = Par1.r[i] - Par2.r[i];
    d = d + c*c;
  }
  return sqrt(d);
}


// Computes the L2 distance between two sets of parameters, including the root.
double parameters_distance(Parameters &Par1, Parameters &Par2) {
  long i,j,k;
  double d, e;
  if (Par1.nedges != Par2.nedges || Par1.nalpha != Par2.nalpha) {
    std::cout << "Error: Can't compare the two parameters." << std::endl;
    exit(1);
  }
 
  d = 0;
  for (i=0; i < Par1.nedges; i++) {
    for(j=0; j < Par1.nalpha; j++) {
      for(k=0; k < Par1.nalpha; k++) {
        e = Par1.tm[i][j][k] - Par2.tm[i][j][k];
        d = d + e*e;
      }
    }
  }

  for (i=0; i < Par1.nalpha; i++) {
    e = Par1.r[i] - Par2.r[i];
    d = d + e*e;
  }

  return sqrt(d);
}



// Prints out the parameters
void print_parameters(Parameters &Par, std::ostream &s) {
  int i,j,k;
  s << Par.nalpha << " " << Par.nedges << std::endl << std::endl;
  for (i=0; i < Par.nalpha; i++) {
    s << Par.r[i] << " ";
  }
  s << std::endl << std::endl;

  for (k=0; k < Par.nedges; k++) {
    for(i=0; i < Par.nalpha; i++) {
      for(j=0; j < Par.nalpha; j++) {
	s << Par.tm[k][i][j] << " ";
      }
      s << std::endl;
    }
    s << std::endl;
  }
  s << std::endl;
}

// Saves the parameters into a file.
void save_parameters(Parameters &Par, std::string const &fname) {
  std::fstream st;
  st.precision(15);

  st.open(fname.c_str(), std::ios::out);
  if (!st.is_open()) {
    std::cout << "Could not open file: " << fname << std::endl;
  }
  print_parameters(Par, st);
  st.close();
}

// Reads the parameters from a file.
void read_parameters(Parameters &Par, std::string const &fname) {
  long nalpha, nedges;
  long i,j,k;

  std::fstream fpar;
  fpar.open(fname.c_str(), std::ios::in);
  fpar >> nalpha >> nedges;
  if(Par.nalpha != nalpha || Par.nedges != nedges) {
    std::cout << "The number of edges or nalpha in the file are not what was expected" << std::endl;
    exit(-1);
  }

  for(i=0; i < nalpha; i++) {
    fpar >> Par.r[i];
  }

  for(i=0; i < nedges; i++) {
    for(j=0; j < nalpha; j++) {
      for(k=0; k < nalpha; k++) {
        fpar >> Par.tm[i][j][k];
      }
    }
  }
  fpar.close();
}


void copy_parameters(Parameters &source, Parameters &target) {
  target.nedges = source.nedges;
  target.nalpha = source.nalpha;
 
  target.tm.resize(source.nedges);
  for (long i=0; i < source.nedges; i++) {
    target.tm[i].resize(source.nalpha);
    for (long j=0; j < source.nalpha; j++) {
      target.tm[i][j].resize(source.nalpha);
      for (long k=0; k < source.nalpha; k++) {
        target.tm[i][j][k] = source.tm[i][j][k];
      }
    }
  }

  target.r.resize(source.nalpha);
  for (long i=0; i < source.nalpha; i++) {
    target.r[i] = source.r[i];
  }
}
