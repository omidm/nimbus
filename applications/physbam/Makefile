
.PHONY: default dependencies_t release debug

default: release

dependencies_t:
	sudo ./install-dependencies.sh

release:
	cd physbam-lib/Public_Library; mkdir -p Build/Release;
	cd physbam-lib/Public_Library/Build/Release; \
    cmake -D PHYSBAM_USE_OPENMPI=ON -D PHYSBAM_USE_LAMMPI=OFF ../..; cmake -LH ../..
	$(MAKE) -C physbam-lib/Public_Library/Build/Release
	cd water/; mkdir -p Build/Release;
	cd water/Build/Release; \
    cmake -D PHYSBAM_USE_OPENMPI=ON -D PHYSBAM_USE_LAMMPI=OFF ../..; cmake -LH ../..
	$(MAKE) -C water/Build/Release
	cd smoke/; mkdir -p Build/Release;
	cd smoke/Build/Release; \
    cmake -D PHYSBAM_USE_OPENMPI=ON -D PHYSBAM_USE_LAMMPI=OFF ../..; cmake -LH ../..
	$(MAKE) -C smoke/Build/Release
	cd physbam-app/opengl_3d/; mkdir -p Build/Release;
	cd physbam-app/opengl_3d/Build/Release; \
    cmake -D PHYSBAM_USE_OPENMPI=ON -D PHYSBAM_USE_LAMMPI=OFF ../..; cmake -LH ../..
	$(MAKE) -C physbam-app/opengl_3d/Build/Release
	cd physbam-app/water/; mkdir -p Build/Release;
	cd physbam-app/water/Build/Release; \
    cmake -D PHYSBAM_USE_OPENMPI=ON -D PHYSBAM_USE_LAMMPI=OFF ../..; cmake -LH ../..
	$(MAKE) -C physbam-app/water/Build/Release

debug:
	cd physbam-lib/Public_Library; mkdir -p Build/Debug;
	cd physbam-lib/Public_Library/Build/Debug; \
    cmake -D PHYSBAM_USE_OPENMPI=ON -D PHYSBAM_USE_LAMMPI=OFF ../..; cmake -LH ../..
	$(MAKE) -C physbam-lib/Public_Library/Build/Debug
	cd water/; mkdir -p Build/Debug;
	cd water/Build/Debug; \
    cmake -D PHYSBAM_USE_OPENMPI=ON -D PHYSBAM_USE_LAMMPI=OFF ../..; cmake -LH ../..
	$(MAKE) -C water/Build/Debug
	cd smoke/; mkdir -p Build/Debug;
	cd smoke/Build/Debug; \
    cmake -D PHYSBAM_USE_OPENMPI=ON -D PHYSBAM_USE_LAMMPI=OFF ../..; cmake -LH ../..
	$(MAKE) -C smoke/Build/Debug
	cd physbam-app/opengl_3d/; mkdir -p Build/Debug;
	cd physbam-app/opengl_3d/Build/Debug; \
    cmake -D PHYSBAM_USE_OPENMPI=ON -D PHYSBAM_USE_LAMMPI=OFF ../..; cmake -LH ../..
	$(MAKE) -C physbam-app/opengl_3d/Build/Debug
	cd physbam-app/water/; mkdir -p Build/Debug;
	cd physbam-app/water/Build/Debug; \
    cmake -D PHYSBAM_USE_OPENMPI=ON -D PHYSBAM_USE_LAMMPI=OFF ../..; cmake -LH ../..
	$(MAKE) -C physbam-app/water/Build/Debug


CLEANSUBDIRS  = physbam-lib/Public_Library/Build/Release
CLEANSUBDIRS += physbam-lib/Public_Library/Build/Debug
CLEANSUBDIRS += physbam-app/water/Build/Release
CLEANSUBDIRS += physbam-app/water/Build/Debug
CLEANSUBDIRS += physbam-app/opengl_3d/Build/Release
CLEANSUBDIRS += physbam-app/opengl_3d/Build/Debug
CLEANSUBDIRS += water/Build/Release
CLEANSUBDIRS += water/Build/Debug

CMAKECACHE  = physbam-lib/External_Libraries/src/fltk-build/CMakeCache.txt
CMAKECACHE += physbam-lib/External_Libraries/src/zlib-build/CMakeCache.txt 
CMAKECACHE += physbam-lib/External_Libraries/src/freeglut-build/CMakeCache.txt
CMAKECACHE += physbam-lib/External_Libraries/src/libjpeg-build/CMakeCache.txt
CMAKECACHE += physbam-lib/External_Libraries/src/libpng-build/CMakeCache.txt

.PHONY: $(CLEANSUBDIRS) $(CMAKECACHE) clean

clean: $(CLEANSUBDIRS) $(CMAKECACHE)

$(CLEANSUBDIRS):
	@if [ -d $@ ]; then \
    $(MAKE) -C $@ clean; \
	  cd $@; cd ../..; rm -rf Build; \
  fi

$(CMAKECACHE):
	@rm -f $@;\

