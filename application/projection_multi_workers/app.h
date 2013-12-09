/*
 * The application specification of PhysBAM projection.
 *
 * Author: Hang Qu <quhang@stanford.edu>
 */

#ifndef __PROJECTION_MAIN__ 
#define __PROJECTION_MAIN__

#include "shared/nimbus.h"

// TODO: Use clause is not allowed in header file.
using nimbus::Data;
using nimbus::Application;

namespace PhysBAM {
template<class TV> class PROJECTION_DRIVER;
template<class T, int d> class VECTOR;
class MPI_WORLD;
} // namespace PhysBAM

class App : public Application {
public:
	App(int rankID) {
		_rankID = rankID;
	}
	void Load();
	// Initializes projection driver.
	void InitMain(int argc, char* argv[]);
	// Tears down the projection driver.
	void FinishMain();
	// Driver contains all the meta-data and functions to do projection on
	// a grid.
	PhysBAM::PROJECTION_DRIVER< PhysBAM::VECTOR<float,2> >* app_driver;
	
	int _rankID;
};

// "ProfileData" contains nothing for now.
class ProfileData : public Data {
public:
	ProfileData() {
	}
	virtual void Create() {
	}
	virtual Data * Clone() {
		return new ProfileData;
	}
};

#endif
