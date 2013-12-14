
#################################
# ENTER application target here #
#################################

TARGET_WATER = worker

##########
# common #
##########

all: executable
	 @echo "Making ${TARGET_WATER}"

NIMBUS_ROOT = ../..
include ${NIMBUS_ROOT}/Makeinclude

H_FILES       = $(wildcard *.h)
CPP_FILES     = $(wildcard *.cc)
OBJECT_FILES  = $(subst .cc,.o,$(CPP_FILES))

PHYSBAM_DIR    = $(NIMBUS_ROOT)/application/physbam-lib/Public_Library/
APP_LIB        = -lwater_app
debug: APP_LIB = -lwater_app-debug-debug

####################################
# ENTER application directory here #
####################################

APP_DIR        = $(NIMBUS_ROOT)/application/water_alternate_coarse/Build/Release
debug: APP_DIR = $(NIMBUS_ROOT)/application/water_alternate_coarse/Build/Debug

####################
# ENTER flags here #
####################

IFLAGS  += -I$(LDIR)
IFLAGS  += -I$(APP_DIR) -I$(PHYSBAM_DIR)
LFLAGS  += -lnimbus -lprotobuf -lboost_thread-mt -lboost_system-mt -lpthread $(APP_LIB)
LDFLAGS += -L$(APP_DIR) -Wl,-rpath $(APP_DIR)

#############################
# add the dependencies here #
#############################

%.o: %.cc
	$(CPP) $(CFLAGS) $(IFLAGS) -c $< -o $@

####################################
# ENTER rules to build target here #
####################################

executable: $(CPP_FILES) $(H_FILES) $(OBJECT_FILES)
	$(CPP) $(CFLAGS) $(IFLAGS) $(LDFLAGS) $(LFLAGS) $(OBJECT_FILES) -o $(TARGET_WATER) $(LFLAGS)

release: executable

debug: executable

clean: clean-files
	\rm -f $(TARGET_WATER)
