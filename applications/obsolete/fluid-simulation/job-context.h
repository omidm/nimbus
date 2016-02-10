#ifndef JOB_CONTEXT
#define JOB_CONTEXT
// Between each step, there are some common context varialbes.
struct JobContext {
  // Temporary: -1 means invalid.
  int sub_step;
  float dt; 
  float target_time;
  bool done_after_this_substep;

  struct JobContext *clone() const {
    struct JobContext *result = new struct JobContext;
    std::memcpy(result, this, sizeof(this));
    return result;
  }
};
#endif  // JOB_CONTEXT

