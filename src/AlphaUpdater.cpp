/////////////////////////////////////////////////////////////
///// Class that regulates rank-one updates to small portions of
///// the giant alpha symmetric matrix, allowing finer-grained locking
/////////////////////////////////////////////////////////////
#include "AlphaUpdater.h"
#include "Std.h"
#include "UseTMV.h"

astrometry::AlphaUpdater::AlphaUpdater(tmv::SymMatrix<double>& alpha_,
				       int nMaps_,
				       int nLocks): alpha(alpha_), nMaps(nMaps_) {
#ifdef _OPENMP
  // Decide how many map elements per lock
  int nBlocks = nMaps*(nMaps+1)/2;
  if (nLocks > nBlocks || nLocks<=0) {
    blockLength = 1;
    nLocks = nBlocks;
  } else {
    blockLength = nBlocks / nLocks;
    nLocks = (nBlocks+blockLength-1) / blockLength;
  }
  locks.resize(nLocks);
  // build array of locks
  for (int i=0; i<locks.size(); i++)
    omp_init_lock(&locks[i]);
#endif
}
astrometry::AlphaUpdater::~AlphaUpdater() {
#ifdef _OPENMP    // Release locks
  for (int i=0; i<locks.size(); i++)
    omp_destroy_lock(&locks[i]);
#endif
}

// Update submatrix with corner at (startIndex1,startIndex2) with += scalar * v1 ^ v2
// map[12] are sequence numbers of map components in the PixelMapCollection, used
// to divide the matrix into blocks for locking in multithreaded case
void astrometry::AlphaUpdater::rankOneUpdate(int map1, int startIndex1,
					     const tmv::ConstVectorView<double>& v1, 
					     int map2, int startIndex2,
					     const tmv::ConstVectorView<double>& v2,
					     double scalar) {
  if (v1.size() <= 0 || v2.size() <=0) return;	// Nothing to add.
  bool diagonal = (map1 == map2);
  bool swapVectors = (map1 < map2);
#ifdef _OPENMP
  // Get the block index corresponding to these 2 keys
  int block = (swapVectors ? (map2*(map2+1)/2 + map1) : (map1*(map1+1)/2 + map2))
    / blockLength;
  Assert(block < locks.size());
  // Set lock for its region - will wait here if busy!
  omp_set_lock(&locks[block]);
  // Pre-multiply the vectors?
#endif
  if (diagonal) {
    tmv::SymMatrix<double> tmp = v1 ^ v1;
    tmp *= scalar;
    alpha.subSymMatrix(startIndex1, startIndex1 + v1.size()) += tmp;
  } else if (swapVectors) {
    alpha.subMatrix(startIndex2, startIndex2+v2.size(),
		    startIndex1, startIndex1+v1.size()) += scalar * v2 ^ v1;
  } else {
    alpha.subMatrix(startIndex1, startIndex1+v1.size(),
		    startIndex2, startIndex2+v2.size()) += scalar * v1 ^ v2;
  }
#ifdef _OPENMP
  omp_unset_lock(&locks[block]);
#endif
}
