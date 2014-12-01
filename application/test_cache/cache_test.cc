#include <ctime>
#include "application/test_cache/cache_app_var.h"
#include "application/test_cache/cache_face_array.h"
#include "application/test_cache/cache_scalar_array.h"
#include "application/test_cache/cache_test.h"
#include "application/test_cache/data_face_array.h"
#include "application/test_cache/data_scalar_array.h"
#include "data/physbam/physbam_data.h"
#include "shared/nimbus.h"

namespace application {

// typedefs
typedef std::vector<nimbus::Data*> DataArray;
typedef typename PhysBAM::FACE_INDEX<3> FaceIndex;
typedef typename PhysBAM::ARRAY<float, FaceIndex> PhysBAMFaceArray;

typedef typename PhysBAM::VECTOR<int, 3> ScalarIndex;
typedef typename PhysBAM::ARRAY<float, ScalarIndex> PhysBAMScalarArray;

// constants
const int scale = 64;
const nimbus::GeometricRegion test_region(1, 1, 1, scale, scale, scale);
const int ghost_width = 3;
const int data_partitions = 3;
const int ghost_start[] = {1, 1 + ghost_width, scale - ghost_width + 1};
const int ghost_ls[] = {ghost_width, scale - 2 * ghost_width, ghost_width};

CacheFaceArray<float> cache_fa_proto(test_region, 0, true, "face_array");
DataFaceArray<float> data_fa_proto("face_array");

CacheScalarArray<float> cache_sa_proto(test_region, 0, true, "scalar_array");
DataScalarArray<float> data_sa_proto("scalar_array");

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
    for (int dim = 1; dim <= 3; ++dim) {
      int scalei = dim==1? scale + 1 : scale;
      int scalej = dim==2? scale + 1 : scale;
      int scalek = dim==3? scale + 1 : scale;
      for (int i = 1; i <= scalei; ++i) {
        for (int j = 1; j <= scalej; ++j) {
          for (int k = 1; k <= scalek; ++k) {
            typename PhysBAM::VECTOR<int, 3> index(i, j, k);
            (*pfa)(dim, index) = dim * (i + j + k);
          }
        }
      }
    }
  }

  if (dynamic_cast< CacheScalarArray<float> * >(cv)) {
    printf("Typecast to scalar array ...\n");
    CacheScalarArray<float> *cf = dynamic_cast< CacheScalarArray<float> * >(cv);
    PhysBAMScalarArray *psa = cf->data();
    for (int i = 1; i <= scale; ++i) {
      for (int j = 1; j <= scale; ++j) {
	for (int k = 1; k <= scale; ++k) {
	  typename PhysBAM::VECTOR<int, 3> index(i, j, k);
	  (*psa)(index) = (i + j + k);
	}
      }
    }
  }

  // Test aggregate copy time
  // BEGIN TIMER
  struct timespec t;
  clock_gettime(CLOCK_REALTIME, &t);
  double start_time = t.tv_sec + .000000001 * static_cast<double>(t.tv_nsec);
  ca->Write(data, test_region);
  // END TIMER
  clock_gettime(CLOCK_REALTIME, &t);
  double end_time = t.tv_sec + .000000001 * static_cast<double>(t.tv_nsec);
  double diff_time = end_time - start_time;
  printf("%s, Aggregate, %f s\n", data_proto.name().c_str(), diff_time);

  // Test total of inidividual copy times
  // BEGIN TIMER
  clock_gettime(CLOCK_REALTIME, &t);
  start_time = t.tv_sec + .000000001 * static_cast<double>(t.tv_nsec);
  for (int i = 0; i < data_partitions * data_partitions * data_partitions - 1; ++i) {
    DataArray data_single;
    data_single.push_back(data[i]);
    ca->Write(data_single, data[i]->region());
  }
  // END TIMER
  clock_gettime(CLOCK_REALTIME, &t);
  end_time = t.tv_sec + .000000001 * static_cast<double>(t.tv_nsec);
  diff_time = end_time - start_time;
  printf("%s, Individual, %f s\n", data_proto.name().c_str(), diff_time);
}

void CacheTest::GetWriteTimes() {
  printf("***********Get Write Times: FACE_ARRAY...***********\n");
  GetWriteTime(cache_fa_proto, data_fa_proto);
  printf("***********Get Write Times: SCALAR_ARRAY...***********\n");
  GetWriteTime(cache_sa_proto, data_sa_proto);
}

} // namespace application
