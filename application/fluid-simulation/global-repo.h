#ifndef GLOBOL_REPO_H_
#define GLOBOL_REPO_H_
namespace PhysBAM {
class WATER_EXAMPLE;
class WATER_DRIVER;

// Given that the control flow is not normal, to store global structure
// here makes life easier.
struct GlobalRepo {
  WATER_EXAMPLE *water_example;
  WATER_DRIVER *water_driver;
};
extern struct GlobalRepo *g_global_repo;
}  // namespace PhysBAM
#endif  // GLOBAL_REPO_H_
