#include "shared/nimbus.h"
#include "projection_app.h"

using nimbus::Data;
using nimbus::Job;
using nimbus::Application;

ProjectionApp::ProjectionApp() {
};

void ProjectionApp::Load() {
    printf("Worker beginning to load application\n");

    /* Declare and initialize data, jobs and policies. */

    RegisterJob("main", new Main(this));
    
    
    printf("Finished creating job and data definitions\n");
    printf("Finished loading application\n");
}

Main::Main(Application *app) {
    set_application(app);
};

Job* Main::Clone() {
    printf("Cloning main job\n");
    return new Main(application());
};

void Main::Execute(std::string params, const DataArray& da)
{
    printf("Begin main\n");
    
    printf("Completed main\n");
};

