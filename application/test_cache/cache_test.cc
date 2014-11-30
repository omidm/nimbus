#include "application/test_cache/cache_face_array.h"
#include "application/test_cache/cache_test.h"
#include "application/test_cache/data_face_array.h"
#include "data/cache/cache_defs.h"
#include "data/cache/cache_object.h"
#include "data/cache/cache_table.h"
#include "data/physbam/physbam_data.h"
#include "shared/geometric_region.h"
#include "shared/nimbus.h"
#include "worker/data.h"

namespace application {

// typedefs
typedef std::vector<nimbus::Data*> DataArray;
typedef typename PhysBAM::FACE_INDEX<3> FaceIndex;
typedef typename PhysBAM::ARRAY<float, FaceIndex> PhysBAMFaceArray;

// constants
const int scale = 32;
const nimbus::GeometricRegion test_region(1, 1, 1, scale, scale, scale);
const int ghost_width = 3;
const int data_partitions = 3;
const int ghost_start[] = {1, 1 + ghost_width, scale - ghost_width + 1};

// face array
void CacheTest::TestCacheFaceArray() {
  printf("Running test cache face array...\n");

  // cache variable
  CacheFaceArray<float> fa_proto(test_region, ghost_width, true, "face_array");
  nimbus::CacheVar* cv = fa_proto.CreateNew(test_region);
  // nimbus objects
  DataArray data;
  for (int i = 0; i < data_partitions; ++i) {
    for (int j = 0; j < data_partitions; ++j) {
      for (int k = 0; k < data_partitions; ++k) {
	// skip central region
	if (i == 1 && j == 1 && k == 1)
	  continue;
        nimbus::Data *d = new DataFaceArray<float>("face_array");
	d->set_region(nimbus::GeometricRegion(ghost_start[i], ghost_start[j], ghost_start[k],
	                                      ghost_width, ghost_width, ghost_width));
        d->Create();
	data.push_back(d);
      }
    }
  }
  CacheFaceArray<float> *cf = dynamic_cast< CacheFaceArray<float> * >(cv);
  // Initialize face array here to something that is not constant.
  PhysBAMFaceArray *pfa = cf->data();
  for (int i = 1; i <= scale; ++i) {
    for (int j = 1; j <= scale; ++j) {
      for (int k = 1; k <= scale; ++k) {
	typename PhysBAM::VECTOR<int, 3> index(i, j, k);
	(*pfa)(1, index) = (i + j + k);
	(*pfa)(2, index) = 2 * (i + j + k);
	(*pfa)(3, index) = 3 * (i + j + k);
      }
    }
  }

  // Test aggregate copy time
  // TODO: BEGIN TIMER
  cf->Write(data, test_region);
  // TODO: END TIMER
  printf("FaceArray, Aggregate, %i\n", 0);

  // Test total of inidividual copy times
  // TODO: BEGIN TIMER
  for (int i = 0; i < data_partitions * data_partitions * data_partitions - 1; ++i) {
    DataArray data_single;
    data_single.push_back(data[i]);
    cf->Write(data_single, data[i]->region());
  }
  // TODO: END TIMER
  printf("FaceArray, Individual, %i\n", 0);

  printf("Exiting from test cache face array...\n");
}

} // namespace application
