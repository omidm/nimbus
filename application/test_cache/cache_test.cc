#include "application/test_cache/cache_app_var.h"
#include "application/test_cache/cache_face_array.h"
#include "application/test_cache/cache_test.h"
#include "application/test_cache/data_face_array.h"
#include "data/physbam/physbam_data.h"
#include "shared/nimbus.h"

namespace application {

// typedefs
typedef std::vector<nimbus::Data*> DataArray;
typedef typename PhysBAM::FACE_INDEX<3> FaceIndex;
typedef typename PhysBAM::ARRAY<float, FaceIndex> PhysBAMFaceArray;

// constants
const int scale = 64;
const nimbus::GeometricRegion test_region(1, 1, 1, scale, scale, scale);
const int ghost_width = 3;
const int data_partitions = 3;
const int ghost_start[] = {1, 1 + ghost_width, scale - ghost_width + 1};
const int ghost_ls[] = {ghost_width, scale - 2 * ghost_width, ghost_width};

CacheFaceArray<float> cache_fa_proto(test_region, 0, true, "face_array");
DataFaceArray<float> data_fa_proto("face_array");

void CacheTest::GetWriteTime(CacheAppVar &cache_proto, nimbus::PhysBAMData &data_proto) {
  printf("Running %s ...\n", data_proto.name().c_str());

  // cache variable
  nimbus::CacheVar* cv = cache_proto.CreateNew(test_region);
  // nimbus objects
  DataArray data;
  for (int i = 0; i < data_partitions; ++i) {
    for (int j = 0; j < data_partitions; ++j) {
      for (int k = 0; k < data_partitions; ++k) {
	// skip central region
	if (i == 1 && j == 1 && k == 1)
	  continue;
        nimbus::Data *d = data_proto.Clone();
	d->set_region(nimbus::GeometricRegion(ghost_start[i], ghost_start[j], ghost_start[k],
	                                      ghost_ls[i], ghost_ls[j], ghost_ls[k]));
        d->Create();
	data.push_back(d);
      }
    }
  }
  CacheAppVar *ca = dynamic_cast< CacheAppVar * >(cv);

  // Initialize physbam structures in the cachge objects to something that
  // is not constant.
  if (dynamic_cast< CacheFaceArray<float> * >(cv)) {
    printf("Typecast to face array ...\n");
    CacheFaceArray<float> *cf = dynamic_cast< CacheFaceArray<float> * >(cv);
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
  }

  // Test aggregate copy time
  // TODO: BEGIN TIMER
  ca->Write(data, test_region);
  // TODO: END TIMER
  printf("%s, Aggregate, %i\n", data_proto.name().c_str(), 0);

  // Test total of inidividual copy times
  // TODO: BEGIN TIMER
  for (int i = 0; i < data_partitions * data_partitions * data_partitions - 1; ++i) {
    DataArray data_single;
    data_single.push_back(data[i]);
    ca->Write(data_single, data[i]->region());
  }
  // TODO: END TIMER
  printf("%s, Individual, %i\n", data_proto.name().c_str(), 0);
}

void CacheTest::GetWriteTimes() {
  GetWriteTime(cache_fa_proto, data_fa_proto);
}

} // namespace application
