/*
 * The job specification of PhysBAM projection for inner loop.
 * MPI is heavily used in this version.
 *
 */

#include <PhysBAM_Tools/Parallel_Computation/MPI_WORLD.h>
#include "projection_driver.h"
#include "projection_example.h"

#include "app.h"
#include "data_impl.h"
#include "job_impl.h"

// TODO The idset implementation should offer the helper function.
#define READ_0() read.clear()
#define READ_1(x) read.clear(); read.insert(x)
#define READ_2(x,y) READ_1(x); read.insert(y)
#define READ_3(x,y,z) READ_2(x,y); read.insert(z)
#define READ_4(x,y,z,e) READ_3(x,y,z); read.insert(e)
#define READ_5(x,y,z,e,f) READ_4(x,y,z,e); read.insert(f)
#define WRITE_0() write.clear()
#define WRITE_1(x) write.clear(); write.insert(x)
#define WRITE_2(x,y) WRITE_1(x); write.insert(y)
#define WRITE_3(x,y,z) WRITE_2(x,y); write.insert(z)
#define BEFORE_0() before.clear();
#define BEFORE_1(x) before.clear(); before.insert(x)
#define BEFORE_2(x,y) BEFORE_1(x); before.insert(y)
#define BEFORE_3(x,y,z) BEFORE_2(x,y); before.insert(z)
#define AFTER_0() after.clear();
#define AFTER_1(x) after.clear(); after.insert(x)
#define AFTER_2(x,y) AFTER_1(x); after.insert(y)
#define AFTER_3(x,y,z) AFTER_2(x,y); after.insert(z)

OneIteration::OneIteration(Application *app) {
	set_application(app);
}

void OneIteration::Execute(Parameter params, const DataArray& da) {
	App* projection_app = dynamic_cast<App*>(application());
	dbg(DBG_PROJ, "||Iteration job starts on worker %d.\n",
			projection_app->_rankID);
	PartialNorm *partial_norm = reinterpret_cast<PartialNorm*>(da[1]);
	PhysBAM::PROJECTION_DRIVER< PhysBAM::VECTOR<float,2> >* app_driver =
			projection_app->app_driver;
	app_driver->projection_internal_data->iteration++;
	app_driver->pcg_mpi->DoPrecondition(app_driver->projection_internal_data,
			app_driver->projection_data);
	app_driver->pcg_mpi->CalculateBeta(app_driver->projection_internal_data,
			app_driver->projection_data);
	app_driver->pcg_mpi->UpdateSearchVector(
			app_driver->projection_internal_data, app_driver->projection_data);
	app_driver->pcg_mpi->ExchangeSearchVector(
			app_driver->projection_internal_data, app_driver->projection_data);
	app_driver->pcg_mpi->UpdateTempVector(app_driver->projection_internal_data,
			app_driver->projection_data);
	app_driver->pcg_mpi->CalculateAlpha(app_driver->projection_internal_data,
			app_driver->projection_data);
	app_driver->pcg_mpi->UpdateOtherVectors(
			app_driver->projection_internal_data, app_driver->projection_data);
	app_driver->pcg_mpi->CalculateResidual(
			app_driver->projection_internal_data, app_driver->projection_data);
	partial_norm->norm_ = app_driver->projection_internal_data->partial_norm;
	dbg(DBG_PROJ, "||Iteration job finishes on worker %d.\n",
			projection_app->_rankID);
}

Job* OneIteration::Clone() {
	return new OneIteration(application());
}
