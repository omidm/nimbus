#ifndef ADVECTION_VELOCITY_H
namespace ADVECT_VELOCITY_NS {
void* advect_velocity_worker(void *arg);
void* advect_velocity_fetcher(void *arg);
void* advect_velocity_fetcher_network(void *arg);
}
#endif
