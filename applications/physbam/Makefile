
.PHONY: dependencies physbam-lib water-physbam water-nimbus opengl_3d

all: physbam-lib water-physbam water-nimbus opengl_3d

dependencies:
	./install-dependencies.sh

physbam-lib: dependencies
	cd physbam-lib/Public_Library; mkdir -p Build/Release;
	cd physbam-lib/Public_Library/Build/Release; \
    cmake -D PHYSBAM_USE_OPENMPI=ON -D PHYSBAM_USE_LAMMPI=OFF ../..
	$(MAKE) -C physbam-lib/Public_Library/Build/Release

water-physbam: dependencies physbam-lib
	cd physbam-app/water/; mkdir -p Build/Release;
	cd physbam-app/water/Build/Release; cmake ../..
	$(MAKE) -C physbam-app/water/Build/Release

water-nimbus: dependencies physbam-lib water-physbam
	cd water/; mkdir -p Build/Release;
	cd water/Build/Release; cmake ../..
	$(MAKE) -C water/Build/Release

opengl_3d: dependencies physbam-lib water-physbam water-nimbus
	cd physbam-app/opengl_3d/; mkdir -p Build/Release;
	cd physbam-app/opengl_3d/Build/Release; cmake ../..
	$(MAKE) -C physbam-app/opengl_3d/Build/Release


.PHONY: physbam-lib-debug water-physbam-debug water-nimbus-debug

physbam-lib-debug: dependencies
	cd physbam-lib/Public_Library; mkdir -p Build/Debug;
	cd physbam-lib/Public_Library/Build/Debug; \
    cmake -D PHYSBAM_USE_OPENMPI=ON -D PHYSBAM_USE_LAMMPI=OFF ../..
	$(MAKE) -C physbam-lib/Public_Library/Build/Debug

water-physbam-debug: dependencies physbam-lib-debug
	cd physbam-app/water/; mkdir -p Build/Debug;
	cd physbam-app/water/Build/Debug; cmake ../..
	$(MAKE) -C physbam-app/water/Build/Debug

water-nimbus-debug: dependencies physbam-lib-debug
	cd water/; mkdir -p Build/Debug;
	cd water/Build/Debug; cmake ../..
	$(MAKE) -C water/Build/Debug


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

.PHONY: $(CLEANSUBDIRS) clean

clean: $(CLEANSUBDIRS) $(CMAKECACHE)

$(CLEANSUBDIRS):
	@if [ -d $@ ]; then \
    $(MAKE) -C $@ clean; \
	  cd $@; cd ../..; rm -rf Build; \
  fi

$(CMAKECACHE):
	@if [ -f $@ ]; then \
    rm $@;\
  fi

