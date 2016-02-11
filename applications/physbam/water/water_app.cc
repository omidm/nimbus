/* Copyright 2013 Stanford University.
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 vd* are met:
 *
 * - Redistributions of source code must retain the above copyright
 *   notice, this list of conditions and the following disclaimer.
 *
 * - Redistributions in binary form must reproduce the above copyright
 *   notice, this list of conditions and the following disclaimer in the
 *   documentation and/or other materials provided with the
 *   distribution.
 *
 * - Neither the name of the copyright holders nor the names of
 *   its contributors may be used to endorse or promote products derived
 *   from this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
 * FOR A PARTICULAR PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL
 * THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT,
 * INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
 * SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
 * HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,
 * STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED
 * OF THE POSSIBILITY OF SUCH DAMAGE.
 */

/*
 * Author: Chinmayee Shah <chinmayee.shah@stanford.edu>
 */

#include "applications/physbam/water//app_utils.h"
#include "applications/physbam/water//app_data_include.h"
#include "applications/physbam/water//app_data_prototypes.h"
#include "applications/physbam/water//data_include.h"
#include "applications/physbam/water//parameters.h"
#include "applications/physbam/water//job_include.h"
#include "applications/physbam/water//reg_def.h"
#include "applications/physbam/water//water_app.h"
#include "src/data/physbam/translator_physbam.h"
#include "src/data/physbam/translator_physbam_old.h"
#include "src/data/scalar_data.h"
#include "src/data/scratch_data_helper.h"
#include <PhysBAM_Tools/Log/DEBUG_SUBSTEPS.h>
#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Tools/Read_Write/Utilities/FILE_UTILITIES.h>
#include "src/shared/dbg.h"
#include "src/shared/geometric_region.h"
#include "src/shared/log.h"
#include "src/shared/nimbus.h"
#include "src/shared/timer.h"
#include "stdio.h"


#include <boost/program_options.hpp>
extern "C" nimbus::Application * ApplicationBuilder(int argc, char *argv[]) {
  namespace po = boost::program_options;

  uint64_t scale;
  uint64_t part_num_x;
  uint64_t part_num_y;
  uint64_t part_num_z;
  uint64_t projection_part_num_x;
  uint64_t projection_part_num_y;
  uint64_t projection_part_num_z;
  uint64_t last_frame;
  uint64_t max_iterations;
  uint64_t iteration_batch;
  uint64_t smart_projection_level;
  float water_level;

  po::options_description desc("Water Application Options");
  desc.add_options()
    ("help,h", "produce help message")

    // Optinal arguments
    ("scale,s", po::value<uint64_t>(&scale), "scale of the simulation") //NOLINT
    ("pnx", po::value<uint64_t>(&part_num_x), "partition number along x") //NOLINT
    ("pny", po::value<uint64_t>(&part_num_y), "partition number along y") //NOLINT
    ("pnz", po::value<uint64_t>(&part_num_z), "partition number along z") //NOLINT
    ("ppnx", po::value<uint64_t>(&projection_part_num_x), "projection partition number along x") //NOLINT
    ("ppny", po::value<uint64_t>(&projection_part_num_y), "projection partition number along y") //NOLINT
    ("ppnz", po::value<uint64_t>(&projection_part_num_z), "projection partition number along z") //NOLINT
    ("last_frame,e", po::value<uint64_t>(&last_frame), "last frame to compute") //NOLINT
    ("maxi", po::value<uint64_t>(&max_iterations), "maximum projection iterations") //NOLINT
    ("ibatch", po::value<uint64_t>(&iteration_batch), "projection iteration batch") //NOLINT
    ("psl", po::value<uint64_t>(&smart_projection_level), "smart projection level") //NOLINT
    ("wl", po::value<float>(&water_level), "initial water level, float between 0 and 1") //NOLINT
    ("dpb", "deactivate projection bottleneck job")
    ("dgw", "deactivate one global write per frame");

  po::variables_map vm;
  try {
    po::store(po::parse_command_line(argc, argv, desc), vm);
  }
  catch(std::exception& e) { // NOLINT
    std::cerr << "ERROR: " << e.what() << "\n";
    return NULL;
  }

  if (vm.count("help")) {
    std::cout << desc << "\n";
    return NULL;
  }

  try {
    po::notify(vm);
  }
  catch(std::exception& e) { // NOLINT
    std::cerr << "ERROR: " << e.what() << "\n";
    return NULL;
  }


  application::WaterApp *app = new application::WaterApp();

  if (vm.count("scale")) {
    app->set_scale(scale);
  }

  if (vm.count("pnx")) {
    app->set_part_num_x(part_num_x);
  }

  if (vm.count("pny")) {
    app->set_part_num_y(part_num_y);
  }

  if (vm.count("pnz")) {
    app->set_part_num_z(part_num_z);
  }

  if (vm.count("ppnx")) {
    app->set_projection_part_num_x(projection_part_num_x);
  }

  if (vm.count("ppny")) {
    app->set_projection_part_num_y(projection_part_num_y);
  }

  if (vm.count("ppnz")) {
    app->set_projection_part_num_z(projection_part_num_z);
  }

  if (vm.count("last_frame")) {
    app->set_last_frame(last_frame);
  }

  if (vm.count("psl")) {
    app->set_smart_projection(smart_projection_level);
  }

  if (vm.count("wl")) {
    app->set_water_level(water_level);
  }

  if (vm.count("maxi")) {
    app->set_max_iterations(max_iterations);
  }

  if (vm.count("ibatch")) {
    app->set_iteration_batch(iteration_batch);
  }

  if (vm.count("dpb")) {
    app->set_spawn_projection_loop_bottleneck(false);
  }

  if (vm.count("dgw")) {
    app->set_global_write(false);
  }

  return app;
}



namespace application {

nimbus::PartitionHandler ph;

uint64_t kScale;
uint64_t kAppPartNum;
uint64_t kAppPartNumX;
uint64_t kAppPartNumY;
uint64_t kAppPartNumZ;
uint64_t kProjAppPartNum;
uint64_t kProjAppPartNumX;
uint64_t kProjAppPartNumY;
uint64_t kProjAppPartNumZ;
nimbus::GeometricRegion kDefaultRegion;
uint64_t kLastFrame;
uint64_t kMaxIterations;
uint64_t kIterationBatch;
bool kSpawnProjectionLoopBottleneck;
bool kUseGlobalWrite;

    WaterApp::WaterApp() {
      scale_ = DEFAULT_SCALE;
      part_num_x_ = DEFAULT_APP_PART_NUM_X;
      part_num_y_ = DEFAULT_APP_PART_NUM_Y;
      part_num_z_ = DEFAULT_APP_PART_NUM_Z;
      projection_part_num_x_ = DEFAULT_APP_PROJ_PART_NUM_X;
      projection_part_num_y_ = DEFAULT_APP_PROJ_PART_NUM_Y;
      projection_part_num_z_ = DEFAULT_APP_PROJ_PART_NUM_Z;
      last_frame_ = DEFAULT_LAST_FRAME;
      max_iterations_ = DEFAULT_MAX_ITERATIONS;
      iteration_batch_ = DEFAULT_ITERATION_BATCH;
      spawn_projection_loop_bottleneck_ =
        DEFAULT_SPAWN_PROJECTION_LOOP_BOTTLENECK;
      global_write_ = DEFAULT_USE_GLOBAL_WRITE;
      translator_log = NULL;
      smart_projection_ = DEFAULT_PROJECTION_SMART; 
      water_level_ = DEFAULT_WATER_LEVEL; 
    }

    float WaterApp::water_level() {
      return water_level_;
    }

    void WaterApp::set_scale(uint64_t scale) {
      scale_ = scale;
    }

    void WaterApp::set_part_num_x(uint64_t part_num_x) {
      part_num_x_ = part_num_x;
    }

    void WaterApp::set_part_num_y(uint64_t part_num_y) {
      part_num_y_ = part_num_y;
    }

    void WaterApp::set_part_num_z(uint64_t part_num_z) {
      part_num_z_ = part_num_z;
    }

    void WaterApp::set_projection_part_num_x(uint64_t projection_part_num_x) {
      projection_part_num_x_ = projection_part_num_x;
    }

    void WaterApp::set_projection_part_num_y(uint64_t projection_part_num_y) {
      projection_part_num_y_ = projection_part_num_y;
    }

    void WaterApp::set_projection_part_num_z(uint64_t projection_part_num_z) {
      projection_part_num_z_ = projection_part_num_z;
    }

    void WaterApp::set_last_frame(uint64_t last_frame) {
      last_frame_ = last_frame;
    }

    void WaterApp::set_max_iterations(uint64_t max_iterations) {
      max_iterations_ = max_iterations;
    }

    void WaterApp::set_iteration_batch(uint64_t iteration_batch) {
      iteration_batch_ = iteration_batch;
    }

    void WaterApp::set_spawn_projection_loop_bottleneck(bool flag) {
      spawn_projection_loop_bottleneck_ = flag;
    }


    void WaterApp::set_global_write(bool flag) {
      global_write_ = flag;
    }

    void WaterApp::set_smart_projection(uint64_t level) {
      smart_projection_ = level;
    }

    void WaterApp::set_water_level(float water_level) {
      water_level_ = water_level;
    }

    /* Register data and job types and initialize constant quantities used by
     * application jobs. */
    void WaterApp::Load() {
        ScopeTimer::Initialize(true);
        //nimbus::Timer::Initialize();
        PhysBAM::Initialize_Particles();
        PhysBAM::Initialize_Read_Write_General_Structures();

        // initialize application specific parameters and constants
        assert((part_num_x_ % projection_part_num_x_) == 0);
        assert((part_num_y_ % projection_part_num_y_) == 0);
        assert((part_num_z_ % projection_part_num_z_) == 0);
        assert((water_level_ >= 0) && (water_level_ <= 1));
        kScale = scale_;
        kAppPartNumX = part_num_x_;
        kAppPartNumY = part_num_y_;
        kAppPartNumZ = part_num_z_;
        kAppPartNum = part_num_x_ * part_num_y_ * part_num_z_;
        kProjAppPartNumX = projection_part_num_x_;
        kProjAppPartNumY = projection_part_num_y_;
        kProjAppPartNumZ = projection_part_num_z_;
        kProjAppPartNum = projection_part_num_x_ * projection_part_num_y_ * projection_part_num_z_;
        kDefaultRegion.Rebuild(1, 1, 1, scale_, scale_, scale_);
        kLastFrame = last_frame_;
        kMaxIterations = max_iterations_;
        kIterationBatch = iteration_batch_;
        kSpawnProjectionLoopBottleneck = spawn_projection_loop_bottleneck_;
        kUseGlobalWrite = global_write_;

        // Initialize the app data prototypes.
        InitializeAppDataPrototypes(kDefaultRegion);

        // dbg_add_mode(APP_LOG_STR);
        // dbg_add_mode(TRANSLATE_STR);

        dbg(APP_LOG, "Loading water application\n");

        // Region constants
        // Old code, using the code generated by python scripts -omidm
        // InitializeRegions();

        ph.AddPartitions("kRegW3",
                         kScale, kScale, kScale, 3, 3, 3, 1, 1, 1,
                         nimbus::PartitionHandler::CENTRAL);

        ph.AddPartitions("kRegW3",
                         kScale, kScale, kScale, 3, 3, 3, 1, 1, 1,
                         nimbus::PartitionHandler::OUTER);

        ph.AddPartitions("kRegW3",
                         kScale, kScale, kScale, 3, 3, 3, 1, 1, 1,
                         nimbus::PartitionHandler::INNER);

        ph.AddPartitions("kRegW1",
                         kScale, kScale, kScale, 1, 1, 1, 1, 1, 1,
                         nimbus::PartitionHandler::CENTRAL);

        ph.AddPartitions("kRegW1",
                         kScale, kScale, kScale, 1, 1, 1, 1, 1, 1,
                         nimbus::PartitionHandler::OUTER);

        ph.AddPartitions("kRegW1",
                         kScale, kScale, kScale, 1, 1, 1, 1, 1, 1,
                         nimbus::PartitionHandler::INNER);

        ph.AddPartitions("kRegW0",
                         kScale, kScale, kScale, 0, 0, 0, 1, 1, 1,
                         nimbus::PartitionHandler::CENTRAL);

        ph.AddPartitions("kRegY2W3",
                         kScale, kScale, kScale, 3, 3, 3,
                         kAppPartNumX, kAppPartNumY, kAppPartNumZ,
                         nimbus::PartitionHandler::CENTRAL);

        ph.AddPartitions("kRegY2W3",
                         kScale, kScale, kScale, 3, 3, 3,
                         kAppPartNumX, kAppPartNumY, kAppPartNumZ,
                         nimbus::PartitionHandler::OUTER);

        ph.AddPartitions("kRegY2W3",
                         kScale, kScale, kScale, 3, 3, 3,
                         kAppPartNumX, kAppPartNumY, kAppPartNumZ,
                         nimbus::PartitionHandler::INNER);

        ph.AddPartitions("kRegY2W3",
                         kScale, kScale, kScale, 3, 3, 3,
                         kAppPartNumX, kAppPartNumY, kAppPartNumZ,
                         nimbus::PartitionHandler::CENTRAL_WGB);

        ph.AddPartitions("kRegY2W1",
                         kScale, kScale, kScale, 1, 1, 1,
                         kAppPartNumX, kAppPartNumY, kAppPartNumZ,
                         nimbus::PartitionHandler::CENTRAL);

        ph.AddPartitions("kRegY2W1",
                         kScale, kScale, kScale, 1, 1, 1,
                         kAppPartNumX, kAppPartNumY, kAppPartNumZ,
                         nimbus::PartitionHandler::OUTER);

        ph.AddPartitions("kRegY2W1",
                         kScale, kScale, kScale, 1, 1, 1,
                         kAppPartNumX, kAppPartNumY, kAppPartNumZ,
                         nimbus::PartitionHandler::INNER);

        ph.AddPartitions("kRegY2W1",
                         kScale, kScale, kScale, 1, 1, 1,
                         kAppPartNumX, kAppPartNumY, kAppPartNumZ,
                         nimbus::PartitionHandler::CENTRAL_WGB);

        ph.AddPartitions("kRegY2W0",
                         kScale, kScale, kScale, 0, 0, 0,
                         kAppPartNumX, kAppPartNumY, kAppPartNumZ,
                         nimbus::PartitionHandler::CENTRAL);

        ph.AddPartitions("kRegY2W8",
                         kScale, kScale, kScale, 8, 8, 8,
                         kAppPartNumX, kAppPartNumY, kAppPartNumZ,
                         nimbus::PartitionHandler::INNER);

        ph.AddPartitions("kProjRegY2W3",
                         kScale, kScale, kScale, 3, 3, 3,
                         kProjAppPartNumX, kProjAppPartNumY, kProjAppPartNumZ,
                         nimbus::PartitionHandler::CENTRAL);

        ph.AddPartitions("kProjRegY2W3",
                         kScale, kScale, kScale, 3, 3, 3,
                         kProjAppPartNumX, kProjAppPartNumY, kProjAppPartNumZ,
                         nimbus::PartitionHandler::OUTER);

        ph.AddPartitions("kProjRegY2W3",
                         kScale, kScale, kScale, 3, 3, 3,
                         kProjAppPartNumX, kProjAppPartNumY, kProjAppPartNumZ,
                         nimbus::PartitionHandler::INNER);

        ph.AddPartitions("kProjRegY2W3",
                         kScale, kScale, kScale, 3, 3, 3,
                         kProjAppPartNumX, kProjAppPartNumY, kProjAppPartNumZ,
                         nimbus::PartitionHandler::CENTRAL_WGB);

        ph.AddPartitions("kProjRegY2W1",
                         kScale, kScale, kScale, 1, 1, 1,
                         kProjAppPartNumX, kProjAppPartNumY, kProjAppPartNumZ,
                         nimbus::PartitionHandler::CENTRAL);

        ph.AddPartitions("kProjRegY2W1",
                         kScale, kScale, kScale, 1, 1, 1,
                         kProjAppPartNumX, kProjAppPartNumY, kProjAppPartNumZ,
                         nimbus::PartitionHandler::OUTER);

        ph.AddPartitions("kProjRegY2W1",
                         kScale, kScale, kScale, 1, 1, 1,
                         kProjAppPartNumX, kProjAppPartNumY, kProjAppPartNumZ,
                         nimbus::PartitionHandler::INNER);

        ph.AddPartitions("kProjRegY2W1",
                         kScale, kScale, kScale, 1, 1, 1,
                         kProjAppPartNumX, kProjAppPartNumY, kProjAppPartNumZ,
                         nimbus::PartitionHandler::CENTRAL_WGB);

        ph.AddPartitions("kProjRegY2W0",
                         kScale, kScale, kScale, 0, 0, 0,
                         kProjAppPartNumX, kProjAppPartNumY, kProjAppPartNumZ,
                         nimbus::PartitionHandler::CENTRAL);



        // PhysBAM logging and R/W
        PhysBAM::LOG::Initialize_Logging(false, false, 1<<30, true, kThreadsNum);
        PhysBAM::FILE_UTILITIES::Create_Directory(kOutputDir+"/common");
        PhysBAM::LOG::Instance()->Copy_Log_To_File(kOutputDir+"/common/log.txt", false);

        // REGISTER DATA
        {
            dbg(APP_LOG, "Registering %s\n", APP_FACE_VEL);
            RegisterData(APP_FACE_VEL, new DataFaceArray<float>(APP_FACE_VEL));
            dbg(APP_LOG, "Registering %s\n", APP_FACE_VEL_GHOST);
            RegisterData(APP_FACE_VEL_GHOST, new DataFaceArray<float>(APP_FACE_VEL_GHOST));
            dbg(APP_LOG, "Registering %s\n", APP_PHI);
            RegisterData(APP_PHI, new DataScalarArray<float>(APP_PHI));
            dbg(APP_LOG, "Registering %s\n", APP_POS_PARTICLES);
            RegisterData(APP_POS_PARTICLES, new DataParticleArray(APP_POS_PARTICLES));
            dbg(APP_LOG, "Registering %s\n", APP_NEG_PARTICLES);
            RegisterData(APP_NEG_PARTICLES, new DataParticleArray(APP_NEG_PARTICLES));
            dbg(APP_LOG, "Registering %s\n", APP_POS_REM_PARTICLES);
            RegisterData(APP_POS_REM_PARTICLES, new DataParticleArray(APP_POS_REM_PARTICLES));
            dbg(APP_LOG, "Registering %s\n", APP_NEG_REM_PARTICLES);
            RegisterData(APP_NEG_REM_PARTICLES, new DataParticleArray(APP_NEG_REM_PARTICLES));
            dbg(APP_LOG, "Registering %s\n", APP_LAST_UNIQUE_PARTICLE_ID);
            RegisterData(APP_LAST_UNIQUE_PARTICLE_ID, new nimbus::ScalarData<int>(APP_LAST_UNIQUE_PARTICLE_ID));
            dbg(APP_LOG, "Registering %s\n", APP_DT);
            RegisterData(APP_DT, new nimbus::ScalarData<float>(APP_DT));

            // These Nimbus data types are used in projection but not used in
            // internal projection loop. They are generally used to generate the
            // matrixes or vectors used in internal projection loop.
            // PSI_D.
            dbg(APP_LOG, "Registering %s\n", APP_PSI_D);
            RegisterData(APP_PSI_D, new DataScalarArray<bool>(APP_PSI_D));
            // PSI_N.
            dbg(APP_LOG, "Registering %s\n", APP_PSI_N);
            RegisterData(APP_PSI_N, new DataFaceArray<bool>(APP_PSI_N));
            // PRESSURE.
            dbg(APP_LOG, "Registering %s\n", APP_PRESSURE);
            RegisterData(APP_PRESSURE, new DataScalarArray<float>(APP_PRESSURE));
            // FILLED_REGION_COLORS.
            dbg(APP_LOG, "Registering %s\n", APP_FILLED_REGION_COLORS);
            RegisterData(APP_FILLED_REGION_COLORS,
                new DataScalarArray<int>(APP_FILLED_REGION_COLORS));
            // DIVERGENCE.
            dbg(APP_LOG, "Registering %s\n", APP_DIVERGENCE);
            RegisterData(APP_DIVERGENCE, new DataScalarArray<float>(APP_DIVERGENCE));

            // These Nimbus data types are used in internal projection loop. They
            // are derived from boundary conditions and act as the linkage between
            // inside projection and outside projection.
            // MATRIX_A.
            dbg(APP_LOG, "Registering %s\n", APP_MATRIX_A);
            RegisterData(APP_MATRIX_A, new DataSparseMatrix(APP_MATRIX_A));
            // VECTOR_B.
            dbg(APP_LOG, "Registering %s\n", APP_VECTOR_B);
            RegisterData(APP_VECTOR_B, new DataRawVectorNd(APP_VECTOR_B));
            // VECTOR_X.
            // dbg(APP_LOG, "Registering %s\n", APP_VECTOR_X);
            // RegisterData(APP_VECTOR_X, new DataRawVectorNd(APP_VECTOR_X));
            // INDEX_C2M.
            dbg(APP_LOG, "Registering %s\n", APP_INDEX_C2M);
            RegisterData(APP_INDEX_C2M, new DataRawGridArray(APP_INDEX_C2M));
            // INDEX_M2C.
            dbg(APP_LOG, "Registering %s\n", APP_INDEX_M2C);
            RegisterData(APP_INDEX_M2C, new DataRawArrayM2C(APP_INDEX_M2C));
            // PROJECTION_LOCAL_N.
            dbg(APP_LOG, "Registering %s\n", APP_PROJECTION_LOCAL_N);
            RegisterData(APP_PROJECTION_LOCAL_N,
                new nimbus::ScalarData<int>(APP_PROJECTION_LOCAL_N));
            // PROJECTION_INTERIOR_N.
            dbg(APP_LOG, "Registering %s\n", APP_PROJECTION_INTERIOR_N);
            RegisterData(APP_PROJECTION_INTERIOR_N,
                new nimbus::ScalarData<int>(APP_PROJECTION_INTERIOR_N));

            // These Nimbus data types are used only in internal projection loop.
            // They are mostly static configurations set for projeciton.
            // PROJECTION_LOCAL_TOLERANCE.
            dbg(APP_LOG, "Registering %s\n", APP_PROJECTION_LOCAL_TOLERANCE);
            RegisterData(APP_PROJECTION_LOCAL_TOLERANCE,
                new nimbus::ScalarData<float>(APP_PROJECTION_LOCAL_TOLERANCE));
            // PROJECTION_GLOBAL_TOLERANCE.
            dbg(APP_LOG, "Registering %s\n", APP_PROJECTION_GLOBAL_TOLERANCE);
            RegisterData(APP_PROJECTION_GLOBAL_TOLERANCE,
                new nimbus::ScalarData<float>(APP_PROJECTION_GLOBAL_TOLERANCE));
            // PROJECTION_GLOBAL_N.
            dbg(APP_LOG, "Registering %s\n", APP_PROJECTION_GLOBAL_N);
            RegisterData(APP_PROJECTION_GLOBAL_N,
                new nimbus::ScalarData<int>(APP_PROJECTION_GLOBAL_N));
            // PROJECTION_DESIRED_N.
            dbg(APP_LOG, "Registering %s\n", APP_PROJECTION_DESIRED_ITERATIONS);
            RegisterData(APP_PROJECTION_DESIRED_ITERATIONS,
                new nimbus::ScalarData<int>(APP_PROJECTION_DESIRED_ITERATIONS));

            // These Nimbus data types are used in internal projection loop.
            // PROJECTION_LOCAL_RESIDUAL.
            dbg(APP_LOG, "Registering %s\n", APP_PROJECTION_LOCAL_RESIDUAL);
            RegisterData(APP_PROJECTION_LOCAL_RESIDUAL,
                new nimbus::ScalarData<double>(APP_PROJECTION_LOCAL_RESIDUAL));
            // PROJECTION_LOCAL_RHO.
            dbg(APP_LOG, "Registering %s\n", APP_PROJECTION_LOCAL_RHO);
            RegisterData(APP_PROJECTION_LOCAL_RHO,
                new nimbus::ScalarData<double>(APP_PROJECTION_LOCAL_RHO));
            // PROJECTION_GLOBAL_RHO.
            dbg(APP_LOG, "Registering %s\n", APP_PROJECTION_GLOBAL_RHO);
            RegisterData(APP_PROJECTION_GLOBAL_RHO,
                new nimbus::ScalarData<double>(APP_PROJECTION_GLOBAL_RHO));
            // PROJECTION_GLOBAL_RHO_OLD.
            dbg(APP_LOG, "Registering %s\n", APP_PROJECTION_GLOBAL_RHO_OLD);
            RegisterData(APP_PROJECTION_GLOBAL_RHO_OLD,
                new nimbus::ScalarData<double>(APP_PROJECTION_GLOBAL_RHO_OLD));
            // PROJECTION_LOCAL_DOT_PRODUCT_FOR_ALPHA.
            dbg(APP_LOG, "Registering %s\n",
                APP_PROJECTION_LOCAL_DOT_PRODUCT_FOR_ALPHA);
            RegisterData(APP_PROJECTION_LOCAL_DOT_PRODUCT_FOR_ALPHA,
                new nimbus::ScalarData<double>(
                    APP_PROJECTION_LOCAL_DOT_PRODUCT_FOR_ALPHA));
            // PROJECTION_ALPHA.
            dbg(APP_LOG, "Registering %s\n", APP_PROJECTION_ALPHA);
            RegisterData(APP_PROJECTION_ALPHA,
                new nimbus::ScalarData<float>(APP_PROJECTION_ALPHA));
            // PROJECTION_BETA.
            dbg(APP_LOG, "Registering %s\n", APP_PROJECTION_BETA);
            RegisterData(APP_PROJECTION_BETA,
                new nimbus::ScalarData<float>(APP_PROJECTION_BETA));
            // MATRIX_C.
            dbg(APP_LOG, "Registering %s\n", APP_MATRIX_C);
            RegisterData(APP_MATRIX_C, new DataSparseMatrix(APP_MATRIX_C));
            // VECTOR_Z.
            dbg(APP_LOG, "Registering %s\n", APP_VECTOR_PRESSURE);
            RegisterData(APP_VECTOR_PRESSURE,
                         new DataRawVectorNd(APP_VECTOR_PRESSURE));
            // VECTOR_Z.
            dbg(APP_LOG, "Registering %s\n", APP_VECTOR_Z);
            RegisterData(APP_VECTOR_Z, new DataRawVectorNd(APP_VECTOR_Z));
            // VECTOR_P_META_FORMAT.
            dbg(APP_LOG, "Registering %s\n", APP_VECTOR_P_META_FORMAT);
            RegisterData(APP_VECTOR_P_META_FORMAT,
                         new DataCompressedScalarArray<float>(APP_VECTOR_P_META_FORMAT));
            // VECTOR_TEMP.
            dbg(APP_LOG, "Registering %s\n", APP_VECTOR_TEMP);
            RegisterData(APP_VECTOR_TEMP, new DataRawVectorNd(APP_VECTOR_TEMP));


            dbg(APP_LOG, "Registering scratch %s\n", APP_POS_PARTICLES);
            kScratchPosParticles.RegisterScratchNames(this, new DataParticleArray(APP_POS_PARTICLES));
            dbg(APP_LOG, "Registering scratch %s\n", APP_NEG_PARTICLES);
            kScratchNegParticles.RegisterScratchNames(this, new DataParticleArray(APP_NEG_PARTICLES));
            dbg(APP_LOG, "Registering scratch %s\n", APP_POS_REM_PARTICLES);
            kScratchPosRemParticles.RegisterScratchNames(this, new DataParticleArray(APP_POS_REM_PARTICLES));
            dbg(APP_LOG, "Registering scratch %s\n", APP_NEG_REM_PARTICLES);
            kScratchNegRemParticles.RegisterScratchNames(this, new DataParticleArray(APP_NEG_REM_PARTICLES));
        }

        // APP DATA INNER REGIONS
        {
            kAppDataFaceVel.set_inner_delta(ph.map()["kRegY2W3Inner"][0].Delta());
            kAppDataFaceVelGhost.set_inner_delta(ph.map()["kRegY2W3Inner"][0].Delta());
            kAppDataPsiN.set_inner_delta(ph.map()["kRegY2W1Inner"][0].Delta());
            kAppDataPsiD.set_inner_delta(ph.map()["kRegY2W1Inner"][0].Delta());
            kAppDataPhi3.set_inner_delta(ph.map()["kRegY2W8Inner"][0].Delta());
            kAppDataPhi7.set_inner_delta(ph.map()["kRegY2W8Inner"][0].Delta());
            kAppDataPhi8.set_inner_delta(ph.map()["kRegY2W8Inner"][0].Delta());
            kAppDataPLE.set_inner_delta(ph.map()["kRegY2W3Inner"][0].Delta());
            kAppDataPressure.set_inner_delta(ph.map()["kRegY2W1Inner"][0].Delta());
            kAppDataColors.set_inner_delta(ph.map()["kProjRegY2W1Inner"][0].Delta());
            kAppDataDivergence.set_inner_delta(ph.map()["kProjRegY2W1Inner"][0].Delta());
            kAppDataIndexC2M.set_inner_delta(ph.map()["kProjRegY2W0Central"][0].Delta());
            kAppDataArrayM2C.set_inner_delta(ph.map()["kProjRegY2W0Central"][0].Delta());
            kAppDataMetaP.set_inner_delta(ph.map()["kProjRegY2W1Inner"][0].Delta());
            kAppDataSparseMatrixA.set_inner_delta(ph.map()["kProjRegY2W0Central"][0].Delta());
            kAppDataSparseMatrixC.set_inner_delta(ph.map()["kProjRegY2W0Central"][0].Delta());
            kAppDataVectorB.set_inner_delta(ph.map()["kProjRegY2W0Central"][0].Delta());
            kAppDataVectorPressure.set_inner_delta(ph.map()["kProjRegY2W0Central"][0].Delta());
            kAppDataVectorZ.set_inner_delta(ph.map()["kProjRegY2W0Central"][0].Delta());
            kAppDataVectorTemp.set_inner_delta(ph.map()["kProjRegY2W0Central"][0].Delta());
        }

        // REGISTER JOBS
        {
            RegisterJob(MAIN, new JobMain(this));
            RegisterJob(INITIALIZE, new JobInitialize(this));
            RegisterJob(CALCULATE_DT, new JobCalculateDt(this));
            RegisterJob(UPDATE_GHOST_VELOCITIES, new JobAdjustPhiWithObjects(this));
            RegisterJob(EXTRAPOLATE_PHI, new JobExtrapolatePhi(this));
            RegisterJob(ADVECT_PHI, new JobAdvectPhi(this));
            RegisterJob(STEP_PARTICLES, new JobStepParticles(this));
            RegisterJob(ADVECT_REMOVED_PARTICLES, new JobAdvectRemovedParticles(this));
            RegisterJob(ADVECT_V, new JobAdvectV(this));
            RegisterJob(APPLY_FORCES, new JobApplyForces(this));
            RegisterJob(LOOP_ITERATION, new JobLoopIteration(this));
            RegisterJob(LOOP_ITERATION_PART_TWO, new JobLoopIterationPartTwo(this));
            RegisterJob(LOOP_FRAME, new JobLoopFrame(this));
            RegisterJob(WRITE_OUTPUT, new JobWriteOutput(this));
            RegisterJob(RESEED_PARTICLES, new JobReseedParticles(this));
            RegisterJob(MODIFY_LEVELSET_PART_ONE, new JobModifyLevelsetPartOne(this));
            RegisterJob(MODIFY_LEVELSET_PART_TWO, new JobModifyLevelsetPartTwo(this));
            RegisterJob(MAKE_SIGNED_DISTANCE, new JobMakeSignedDistance(this));
            RegisterJob(ADJUST_PHI, new JobAdjustPhi(this));
            RegisterJob(DELETE_PARTICLES, new JobDeleteParticles(this));
            RegisterJob(REINCORPORATE_PARTICLES, new JobReincorporateRemovedParticles(this));
            RegisterJob(PROJECTION_MAIN, new JobProjectionMain(this));
            RegisterJob(PROJECTION_TRANSFORM_PRESSURE,
                        new JobProjectionTransformPressure(this));
            RegisterJob(PROJECTION_CALCULATE_BOUNDARY_CONDITION_PART_ONE,
                        new JobProjectionCalculateBoundaryConditionPartOne(this));
            RegisterJob(PROJECTION_CALCULATE_BOUNDARY_CONDITION_PART_TWO,
                        new JobProjectionCalculateBoundaryConditionPartTwo(this));
            RegisterJob(PROJECTION_CONSTRUCT_MATRIX,
                        new JobProjectionConstructMatrix(this));
            RegisterJob(PROJECTION_WRAPUP, new JobProjectionWrapup(this));
            RegisterJob(EXTRAPOLATION, new JobExtrapolation(this));
            RegisterJob(SYNCHRONIZE_PARTICLES, new JobSynchronizeParticles(this));

            RegisterJob(PROJECTION_GLOBAL_INITIALIZE,
                        new JobProjectionGlobalInitialize(this));
            RegisterJob(PROJECTION_LOCAL_INITIALIZE,
                        new JobProjectionLocalInitialize(this));
            RegisterJob(PROJECTION_LOOP_ITERATION, new JobProjectionLoopIteration(this));
            RegisterJob(PROJECTION_LOOP_BOTTLENECK, new JobProjectionLoopBottleneck(this));
            RegisterJob(PROJECTION_STEP_ONE, new JobProjectionStepOne(this));
            RegisterJob(PROJECTION_REDUCE_RHO, new JobProjectionReduceRho(this));
            RegisterJob(PROJECTION_STEP_TWO, new JobProjectionStepTwo(this));
            RegisterJob(PROJECTION_STEP_THREE, new JobProjectionStepThree(this));
            RegisterJob(PROJECTION_REDUCE_ALPHA,
                        new JobProjectionReduceAlpha(this));
            RegisterJob(PROJECTION_STEP_FOUR, new JobProjectionStepFour(this));
            RegisterJob(BARRIER_JOB, new JobBarrier(this));

            RegisterStaticConfigPrototype(
                STATIC_CONFIG_VALID_MASK, new StaticConfigValidMask(kDefaultRegion));
            RegisterStaticConfigPrototype(
                STATIC_CONFIG_U_INTERFACE, new StaticConfigUInterface(kDefaultRegion));
            RegisterStaticConfigPrototype(
                STATIC_CONFIG_FORCE, new StaticConfigForce(kDefaultRegion));
            RegisterStaticConfigPrototype(
                STATIC_CONFIG_COLLISION_BODY,
                new StaticConfigCollisionBody(kDefaultRegion));
        }

        nimbus::TranslatorPhysBAM<float>::log = translator_log;
        nimbus::TranslatorPhysBAMOld<PhysBAM::VECTOR<float, 3> >::log = translator_log;
        dbg(APP_LOG, "Completed loading water application\n");
    }

} // namespace application
