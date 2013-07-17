#include "application/fluid-simulation/app.h"

#include <cstdio>
#include <string>

namespace {
void GetJobFromScheduler(std::string *job) {
  *job = "Run fluid simulation";
}
} // namespace

int main() {
  std::string job;
  GetJobFromScheduler(&job);
  ExecuteJob(job);
  return 0;
}
