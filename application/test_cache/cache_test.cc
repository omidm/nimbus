#include "application/test_cache/cache_face_array.h"
// #include "data/physbam/translator_physbam_old.h"
#include "data/cache/cache_defs.h"
#include "data/cache/cache_object.h"
#include "data/cache/cache_table.h"
#include "shared/nimbus.h"
#include "shared/geometric_region.h"
#include "data/physbam/physbam_data.h"
#include "worker/data.h"
#include "application/test_cache/cache_test.h"

namespace application {

typedef std::vector<nimbus::Data*> DataArray;

  void CacheTest::TestCacheFaceArray() {
    const int kScale = 256;
    const nimbus::GeometricRegion kDefaultRegion(1, 1, 1, kScale, kScale, kScale);
    const nimbus::GeometricRegion region(0, 0, 0, 66, 66, 66);
    int ghost_width = 3;
    CacheFaceArray<int> prototype(kDefaultRegion, ghost_width, true, "face_array");
    nimbus::CacheVar* cv = prototype.CreateNew(region);
    DataArray write_set;
    nimbus::PhysBAMData d = nimbus::PhysBAMData();
    char *buffer = (char *)malloc(8);
    d.set_buffer(buffer, 8);
    for (int i = 0; i < kScale; ++i) {
      for (int j = 0; j < kScale; ++j) {
	for (int k = 0; k < kScale; ++k) {
	  nimbus::Data *data = d.Clone();
	  write_set.push_back(data);
	}
      }
    }
    cv->WriteFromCache(write_set, region);
    printf("Exiting from test cache face array...\n");
  }

} // namespace application
