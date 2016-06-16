/*
 *  permutation.cpp
 *  
 *  Created by Ania M. Kedzierska on 11/11/11.
 *  Copyright 2011 Politecnic University of Catalonia, Center for Genomic Regulation.  This is program can be redistributed, modified or else as given by the terms of the GNU General Public License. 
 *  
 */


#include "permutation.h"
#include "model.h"

void permute_rows(TMatrix &tm, Permutation &perm) {
  TMatrix tma;
  tma.resize(tm.size());
  for(unsigned long i=0; i < tm.size(); i++) {
    tma[i].resize(tm[i].size());
  }

  for(unsigned long i=0; i < tm.size(); i++) {
    for(unsigned long j=0; j < tm[i].size(); j++) {
      tma[perm[i]][j] = tm[i][j]; 
    }
  }

  for(unsigned long i=0; i < tm.size(); i++) {
    for(unsigned long j=0; j < tm[i].size(); j++) {
      tm[i][j] = tma[i][j]; 
    }
  }
}


void permute_cols(TMatrix &tm, Permutation &perm) {
  TMatrix tma;
  tma.resize(tm.size());
  for(unsigned long i=0; i < tm.size(); i++) {
    tma[i].resize(tm[i].size());
  }

  for(unsigned long i=0; i < tm.size(); i++) {
    for(unsigned long j=0; j < tm[i].size(); j++) {
      tma[i][perm[j]] = tm[i][j]; 
    }
  }

  for(unsigned long i=0; i < tm.size(); i++) {
    for(unsigned long j=0; j < tm[i].size(); j++) {
      tm[i][j] = tma[i][j]; 
    }
  }
}


void permute_root(Root &r, Permutation &perm) {
  Root ra;
  ra.resize(r.size());
  for(unsigned long i=0; i < r.size(); i++) {
    ra[perm[i]] = r[i];
  }
  for(unsigned long i=0; i < r.size(); i++) {
    r[i] = ra[i];
  }
}




void permute_at_node(Tree &T, long n, Parameters &Par, Permutation &perm) {
  for (long e=0; e < T.nedges; e++) {
    // if n is source
    if (T.edges[e].s == n) {
      permute_rows(Par.tm[e], perm);
    }

    // if n is target
    if (T.edges[e].t == n) {
      permute_cols(Par.tm[e], perm);
    }
  }

  // if n is the root
  if (n == T.nleaves) {
    permute_root(Par.r, perm);
  }
}


void guess_permutation_rec(Tree &T, Model &Mod, Parameters &Par, std::list<Permutation> &L, long node) {
  
  std::list<long> Lout;
  std::list<Permutation>::iterator i;
  std::list<long>::iterator j;
  long k;

  double score, scoremax;
  Permutation permmax;

  permmax.resize(T.nalpha);

  Parameters Par2;

  list_outgoing_edges(T, node, Lout);

  // It's a leaf, so there is nothing to do.
  if (Lout.empty()) return;

  // Fix permutation on descendant subtrees first.
  for (j=Lout.begin(); j != Lout.end(); j++) {
    guess_permutation_rec(T, Mod, Par, L, T.edges[*j].t);
  }

  // Now fix the permutation on a node.
  scoremax = 0;
  for (i = L.begin(); i != L.end(); i++) {

    copy_parameters(Par, Par2);
    permute_at_node(T, node, Par2, *i);
    score = 0;
    for (j=Lout.begin(); j != Lout.end(); j++) {
      for(k=0; k < T.nalpha; k++) {
        score = score + Par2.tm[*j][k][k];
      }
    }

    if (score > scoremax) {
      scoremax = score;
      for(k=0; k < T.nalpha; k++) {
        permmax[k] = (*i)[k];
      }
    }
  }
  permute_at_node(T, node, Par, permmax);

}


void guess_permutation(Tree &T, Model &Mod, Parameters &Par) {

  std::list<Permutation> Lperm;
  (*Mod.list_permutations)(Lperm);

  // Start recursive process on the root.
  guess_permutation_rec(T, Mod, Par, Lperm, T.nleaves);
}



bool is_permutation(Permutation &perm) {
  unsigned long i, j;
  for(i=0; i < perm.size(); i++) {
    if (perm[i] < 0 || perm[i] >= (long)perm.size()) return false;
  }

  for(i=0; i < perm.size(); i++) {
    for(j=i; j < perm.size(); j++) {
      if (i != j && perm[i] == perm[j]) return false;
    }
  }
  return true;
}



long permutation_sign(Permutation &perm) {
  unsigned long i;
  long k, trans;
  Permutation perm2;
  perm2.resize(perm.size());

  for(i=0; i < perm.size(); i++) {
    perm2[i] = perm[i];
  }

  trans = 0;
  bool repeat = true;
  while(repeat) {
    repeat = false;
    for(i=0; i < perm2.size(); i++) {
      k = perm2[i];
      if (k < 0 || k >= (long)perm2.size()) return 0;
      if (k != (long)i) {
        if (perm2[k] != k) {
          perm2[i] = perm2[k];
          perm2[k] = k;
          trans++;
          repeat = true;
        } else {    // If perm2[k] = k, this means that perm2 was not a permutation.
          return 0;
        }
      }
    }
  }
  if (trans % 2 == 0) return 1;
  else return -1;
}
