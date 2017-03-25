// A class that permits efficient parallel construction of symmetric matrix as sum of outer product
// when the outer products are nonzero over short ranges of parameters.
// Each short range of parameters belongs to one "map" and there is a locking mechanism to
// avoid having more than one thread trying to update a given submatrix at the same time.
#ifndef ALPHAUPDATER_H
#define ALPHAUPDATER_H

#ifdef _OPENMP
#include <omp.h>
#endif

#include "UseTMV.h"
#include <vector>

namespace astrometry {
  class AlphaUpdater {
  public:
    AlphaUpdater(tmv::SymMatrix<double>& alpha_,
		 int nMaps_,
		 int nLocks=0);
    ~AlphaUpdater();
    // Update submatrix with corner at (startIndex1,startIndex2) with += scalar * v1 ^ v2
    // map[12] are sequence numbers of map components in the PixelMapCollection, used
    // to divide the matrix into blocks for locking in multithreaded case
    void rankOneUpdate(int map1, int startIndex1, const tmv::ConstVectorView<double>& v1, 
		       int map2, int startIndex2, const tmv::ConstVectorView<double>& v2,
		       double scalar=1.);
  private:
    tmv::SymMatrix<double>& alpha;
    int nMaps;
#ifdef _OPENMP
    int nLocks;
    int blockLength;
    // Lock array
    std::vector<omp_lock_t> locks;
#endif
  };
} // end namespace astrometry

#endif  // ALPHAUPDATER_H
