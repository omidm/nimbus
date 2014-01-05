#include "PhysBAM_Tools/Vectors/VECTOR.h"

/*
 * Utilities to print out information for debug.
 * Can work with normal PhysBAM simulations.
 * Author: Hang Qu <quhang@stanford.edu>
 */
#ifndef DEBUG_UTIL_H_
#define DEBUG_UTIL_H_

// If you want to work with Nimbus, please comment out the following two lines.
#define DBG_WARN 1
#define dbg(X, ...) printf(__VA_ARGS__)

namespace PhysBAM {

template <typename T_GRID>
void PrintGrid(const T_GRID& grid) {
  if (grid.dimension == 2) {
    if (grid.Is_MAC_Grid()) {
      dbg(DBG_WARN, "MAC grid{");
    } else {
      dbg(DBG_WARN, "Normal grid{ ");
    }
    dbg(DBG_WARN, "count=(%d,%d);",
           grid.counts(1),
           grid.counts(2));
    dbg(DBG_WARN, "min=(%.1f,%.1f);",
           (float)grid.domain.min_corner(1),
           (float)grid.domain.min_corner(2));
    dbg(DBG_WARN, "max=(%.1f,%.1f) }\n",
           (float)grid.domain.max_corner(1),
           (float)grid.domain.max_corner(2));
  }
  if (grid.dimension == 3) {
    if (grid.Is_MAC_Grid()) {
      dbg(DBG_WARN, "MAC grid{");
    } else {
      dbg(DBG_WARN, "Normal grid{ ");
    }
    dbg(DBG_WARN, "count=(%d,%d,%d);",
           grid.counts(1),
           grid.counts(2),
           grid.counts(3));
    dbg(DBG_WARN, "min=(%.1f,%.1f,%.1f);",
           (float)grid.domain.min_corner(1),
           (float)grid.domain.min_corner(2),
           (float)grid.domain.min_corner(3));
    dbg(DBG_WARN, "max=(%.1f,%.1f,%.1f) }\n",
           (float)grid.domain.max_corner(1),
           (float)grid.domain.max_corner(2),
           (float)grid.domain.max_corner(3));
  }
}

template<typename T_GRID>
void PrintParticles(const PARTICLE_LEVELSET_UNIFORM<T_GRID>& container,
                    bool detail = false) {
  dbg(DBG_WARN, "Dump info of particles:\n");
  PrintGrid(container.levelset.grid);
  PrintArrayProfile(container.negative_particles);
  if (detail) {
    RANGE<typename T_GRID::VECTOR_INT> domain(
        container.levelset.grid.Domain_Indices(
            container.number_of_ghost_cells));
    domain.max_corner += T_GRID::VECTOR_INT::All_Ones_Vector();
    for (typename T_GRID::NODE_ITERATOR
         iterator(container.levelset.grid, domain);
         iterator.Valid();iterator.Next()) {
      typename T_GRID::VECTOR_INT block=iterator.Node_Index();
      PARTICLE_LEVELSET_PARTICLES<typename T_GRID::VECTOR_T>* cell_particles =
          container.positive_particles(block);
      while (cell_particles) {
        for(int k = 1; k <= cell_particles->array_collection->Size(); k++) {
          assert(cell_particles->radius(k) > 0);
          if (T_GRID::VECTOR_INT::dimension == 2) {
            dbg(DBG_WARN, "X(%.1f,%.1f),", cell_particles->X(k)(1),
                   cell_particles->X(k)(2));
          }
          if (T_GRID::VECTOR_INT::dimension == 3) {
            dbg(DBG_WARN, "X(%.1f,%.1f,%.1f),", cell_particles->X(k)(1),
                   cell_particles->X(k)(2), cell_particles->X(k)(3));
          }
        }
        // final_block=levelset.grid.Block_Index(cell_particles->X(k),keep_particles_in_ghost_cells?number_of_ghost_cells+1:1);
        cell_particles = cell_particles->next;
      }
    }
  }
  dbg(DBG_WARN, "\n[End] Dump info of particles:\n");
}

template <typename T, int d>
void PrintArrayProfile(const ARRAY< T,VECTOR<int,d> >& array) {
  if (d == 2) {
    dbg(DBG_WARN, "count=(%d,%d);",
           array.counts(1),
           array.counts(2));
    dbg(DBG_WARN, "min=(%.1f,%.1f);",
           (float)array.domain.min_corner(1),
           (float)array.domain.min_corner(2));
    dbg(DBG_WARN, "max=(%.1f,%.1f) }\n",
           (float)array.domain.max_corner(1),
           (float)array.domain.max_corner(2));
  }
  if (d == 3) {
    dbg(DBG_WARN, "count=(%d,%d,%d);",
           array.counts(1),
           array.counts(2),
           array.counts(3));
    dbg(DBG_WARN, "min=(%.1f,%.1f,%.1f);",
           (float)array.domain.min_corner(1),
           (float)array.domain.min_corner(2),
           (float)array.domain.min_corner(3));
    dbg(DBG_WARN, "max=(%.1f,%.1f,%.1f) }\n",
           (float)array.domain.max_corner(1),
           (float)array.domain.max_corner(2),
           (float)array.domain.max_corner(3));
  }
}

}  // namespace PhysBAM

#endif  // DEBUG_UTIL_H_
