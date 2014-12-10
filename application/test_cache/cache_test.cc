#include <ctime>
#include "application/test_cache/cache_app_var.h"
#include "application/test_cache/cache_face_array.h"
#include "application/test_cache/cache_scalar_array.h"
#include "application/test_cache/cache_test.h"
#include "application/test_cache/data_face_array.h"
#include "application/test_cache/data_scalar_array.h"
#include "data/physbam/physbam_data.h"
#include "shared/nimbus.h"

# define RUN_TEST 100    // number of times to run test
# define SCALE 64        // scale to test
# define GHOST_WIDTH 3   // ghost width

namespace application {

// typedefs
typedef std::vector<nimbus::Data*> DataArray;
typedef typename PhysBAM::FACE_INDEX<3> FaceIndex;
typedef typename PhysBAM::ARRAY<float, FaceIndex> PhysBAMFaceArrayFloat;
typedef typename PhysBAM::ARRAY<bool, FaceIndex> PhysBAMFaceArrayBool;

typedef typename PhysBAM::VECTOR<int, 3> ScalarIndex;
typedef typename PhysBAM::ARRAY<float, ScalarIndex> PhysBAMScalarArrayFloat;
typedef typename PhysBAM::ARRAY<bool, ScalarIndex> PhysBAMScalarArrayBool;

// constants
const int scale = SCALE;
const nimbus::GeometricRegion test_region(1, 1, 1, scale, scale, scale);
const int ghost_width = GHOST_WIDTH;
const int data_partitions = 3;
const int ghost_start[] = {1, 1 + ghost_width, scale - ghost_width + 1};
const int ghost_ls[] = {ghost_width, scale - 2 * ghost_width, ghost_width};

CacheFaceArray<float> cache_fa_float_proto(test_region, 0, true, "face_array_float");
DataFaceArray<float> data_fa_float_proto("face_array_float");
CacheFaceArray<bool> cache_fa_bool_proto(test_region, 0, true, "face_array_bool");
DataFaceArray<bool> data_fa_bool_proto("face_array_bool");

CacheScalarArray<float> cache_sa_float_proto(test_region, 0, true, "scalar_array_float");
DataScalarArray<float> data_sa_float_proto("scalar_array_float");
CacheScalarArray<bool> cache_sa_bool_proto(test_region, 0, true, "scalar_array_bool");
DataScalarArray<bool> data_sa_bool_proto("scalar_array_bool");

void CacheTest::GetReadWriteTime(CacheAppVar &cache_proto, nimbus::PhysBAMData &data_proto) {
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
    CacheFaceArray<float> *cf = dynamic_cast< CacheFaceArray<float> * >(cv);
    PhysBAMFaceArrayFloat *pfa = cf->data();
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
  if (dynamic_cast< CacheFaceArray<bool> * >(cv)) {
    CacheFaceArray<bool> *cf = dynamic_cast< CacheFaceArray<bool> * >(cv);
    PhysBAMFaceArrayBool *pfa = cf->data();
    for (int dim = 1; dim <= 3; ++dim) {
      int scalei = dim==1? scale + 1 : scale;
      int scalej = dim==2? scale + 1 : scale;
      int scalek = dim==3? scale + 1 : scale;
      for (int i = 1; i <= scalei; ++i) {
        for (int j = 1; j <= scalej; ++j) {
          for (int k = 1; k <= scalek; ++k) {
            typename PhysBAM::VECTOR<int, 3> index(i, j, k);
            (*pfa)(dim, index) = (j < scale/2);
          }
        }
      }
    }
  }

  if (dynamic_cast< CacheScalarArray<float> * >(cv)) {
    CacheScalarArray<float> *cf = dynamic_cast< CacheScalarArray<float> * >(cv);
    PhysBAMScalarArrayFloat *psa = cf->data();
    for (int i = 1; i <= scale; ++i) {
      for (int j = 1; j <= scale; ++j) {
	for (int k = 1; k <= scale; ++k) {
	  typename PhysBAM::VECTOR<int, 3> index(i, j, k);
	  (*psa)(index) = (i + j + k);
	}
      }
    }
  }
  if (dynamic_cast< CacheScalarArray<bool> * >(cv)) {
    CacheScalarArray<bool> *cf = dynamic_cast< CacheScalarArray<bool> * >(cv);
    PhysBAMScalarArrayBool *psa = cf->data();
    for (int i = 1; i <= scale; ++i) {
      for (int j = 1; j <= scale; ++j) {
	for (int k = 1; k <= scale; ++k) {
	  typename PhysBAM::VECTOR<int, 3> index(i, j, k);
	  (*psa)(index) = (j < scale/2);
	}
      }
    }
  }

  struct timespec t;
  double start_time, end_time, diff_time;

  // Test aggregate write time
  // BEGIN TIMER
  clock_gettime(CLOCK_REALTIME, &t);
  start_time = t.tv_sec + .000000001 * static_cast<double>(t.tv_nsec);
  for (int run = 0; run < RUN_TEST; ++run)
    ca->Write(data, test_region);
  // END TIMER
  clock_gettime(CLOCK_REALTIME, &t);
  end_time = t.tv_sec + .000000001 * static_cast<double>(t.tv_nsec);
  diff_time = end_time - start_time;
  printf("%s, aggregate write, %f s\n", data_proto.name().c_str(), diff_time);

  // Test total of inidividual write times
  // BEGIN TIMER
  clock_gettime(CLOCK_REALTIME, &t);
  start_time = t.tv_sec + .000000001 * static_cast<double>(t.tv_nsec);
  for (int run = 0; run < RUN_TEST; ++run) {
    for (int i = 0; i < data_partitions * data_partitions * data_partitions - 1; ++i) {
      DataArray data_single;
      data_single.push_back(data[i]);
      ca->Write(data_single, data[i]->region());
    }
  }
  // END TIMER
  clock_gettime(CLOCK_REALTIME, &t);
  end_time = t.tv_sec + .000000001 * static_cast<double>(t.tv_nsec);
  diff_time = end_time - start_time;
  printf("%s, individual write, %f s\n", data_proto.name().c_str(), diff_time);

  // Test aggregate read time
  // BEGIN TIMER
  clock_gettime(CLOCK_REALTIME, &t);
  start_time = t.tv_sec + .000000001 * static_cast<double>(t.tv_nsec);
  for (int run = 0; run < RUN_TEST; ++run)
    ca->Read(data, test_region);
  // END TIMER
  clock_gettime(CLOCK_REALTIME, &t);
  end_time = t.tv_sec + .000000001 * static_cast<double>(t.tv_nsec);
  diff_time = end_time - start_time;
  printf("%s, aggregate read, %f s\n", data_proto.name().c_str(), diff_time);

  // Test total of inidividual read times
  // BEGIN TIMER
  clock_gettime(CLOCK_REALTIME, &t);
  start_time = t.tv_sec + .000000001 * static_cast<double>(t.tv_nsec);
  for (int run = 0; run < RUN_TEST; ++run) {
    for (int i = 0; i < data_partitions * data_partitions * data_partitions - 1; ++i) {
      DataArray data_single;
      data_single.push_back(data[i]);
      ca->Read(data_single, data[i]->region());
    }
  }
  // END TIMER
  clock_gettime(CLOCK_REALTIME, &t);
  end_time = t.tv_sec + .000000001 * static_cast<double>(t.tv_nsec);
  diff_time = end_time - start_time;
  printf("%s, individual read, %f s\n", data_proto.name().c_str(), diff_time);

  printf("\n");
}

void CacheTest::GetReadWriteTimes() {
  printf("*********** Times for FACE_ARRAY_FLOAT ***********\n");
  GetReadWriteTime(cache_fa_float_proto, data_fa_float_proto);
  printf("*********** Times for FACE_ARRAY_BOOL ***********\n");
  GetReadWriteTime(cache_fa_bool_proto, data_fa_bool_proto);
  printf("***********Times for SCALAR_ARRAY _FLOAT ***********\n");
  GetReadWriteTime(cache_sa_float_proto, data_sa_float_proto);
  printf("***********Times for SCALAR_ARRAY _BOOL ***********\n");
  GetReadWriteTime(cache_sa_bool_proto, data_sa_bool_proto);
}

} // namespace application
