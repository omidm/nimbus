#include <pthread.h>
#include "simple_worker.h"
#include "shared/nimbus.h"
#include "shared/nimbus_types.h"
#include "worker/application.h"
#include "application/protobuf_test/ProtoBuf_test.h"

int main(int argc, char *argv[]) {
  port_t listening_port;
  if (argc < 2) {
    std::cout << "ERROR: provide an integer (1 to 4)." <<
      std::endl;
    exit(-1);
  }
  if (*argv[1] == '1') {
    listening_port = WORKER_PORT_1;
  } else if (*argv[1] == '2') {
    listening_port = WORKER_PORT_2;
  } else if (*argv[1] == '3') {
    listening_port = WORKER_PORT_3;
  } else {
    listening_port = WORKER_PORT_4;
  };
  nimbus_initialize();
  printf("Simple Worker is up!");
  TestApp *app = new TestApp();
  SimpleWorker * w = new SimpleWorker(NIMBUS_SCHEDULER_IP, NIMBUS_SCHEDULER_PORT, listening_port, app);
  w->Run();
}


