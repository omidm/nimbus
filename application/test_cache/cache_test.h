#ifndef NIMBUS_APPLICATION_CACHE_TEST_H_
#define NIMBUS_APPLICATION_CACHE_TEST_H_

#include "application/test_cache/cache_app_var.h"
#include "data/physbam/physbam_data.h"

namespace application {

  class CacheTest {
    public:
      static void GetWriteTimes();
    private:
      static void GetWriteTime(CacheAppVar &cache_proto, nimbus::PhysBAMData &data_proto);
  };

} // namespace application    

#endif // NIMBUS_APPLICATION_CACHE_TEST_H_
