#include "job.h"

Job::Job(int id, FuncToRun handler)
{
  this->id = id;
  this->handler = handler;
}
