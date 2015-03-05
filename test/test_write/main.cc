#include <unistd.h>
#include <sys/syscall.h>
#include "shared/fast_log.hh"
void* branch(void* param) {
  nimbus::timer::InitializeTimers();
  for (long i = 0; i < 1e6; ++i) {
    nimbus::timer::StartTimer(nimbus::timer::kTotal);
    nimbus::timer::StopTimer(nimbus::timer::kTotal);
  }
  return NULL;
}
int main() {
  nimbus::timer::InitializeKeys();
  nimbus::timer::InitializeTimers();
  pthread_t handle;
  nimbus::timer::StartTimer(nimbus::timer::kTotal);
  pthread_create(&handle, NULL, branch, NULL);
  pthread_join(handle, NULL);
  nimbus::timer::StopTimer(nimbus::timer::kTotal);
  nimbus::timer::PrintTimerSummary();
  return 0;
}
