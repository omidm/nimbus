#ifndef NIMBUS_TEST_PROTOBUF_TEST_WORKER_H_
#define NIMBUS_TEST_PROTOBUF_TEST_WORKER_H_

// #define DEBUG_MODE

#include <boost/thread.hpp>
#include <string>
#include <vector>
#include <map>
#include "shared/scheduler_client.h"
#include "shared/serialized_data.h"
#include "shared/cluster.h"
#include "worker/data.h"
#include "worker/job.h"
#include "worker/application.h"
#include "shared/parser.h"
#include "shared/log.h"
#include "worker/worker.h"


class SimpleWorker : public Worker {
  public:
    SimpleWorker(std::string scheduler_ip, port_t scheduler_port,
        port_t listening_port, Application * a);
};





#endif  // NIMBUS_TEST_PROTOBUF_TEST_WORKER_H_
