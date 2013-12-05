#include "./simple_worker.h"

SimpleWorker::SimpleWorker(std::string scheduler_ip, port_t scheduler_port,
        port_t listening_port, Application * a)
: Worker(scheduler_ip, scheduler_port, listening_port, a) {
}

