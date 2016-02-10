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
  // Untyped to be able to avoid unnecessary dependency.
  // [TODO] This should be handled by data map, but be here first for testing.
  void *face_velocities_ghost;
};
extern struct GlobalRepo *g_global_repo;
}  // namespace PhysBAM
#endif  // GLOBAL_REPO_H_
