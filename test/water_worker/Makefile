
#################################
# ENTER application target here #
#################################

TARGET_WATER_SINGLE = worker_single
TARGET_WATER_MULTIPLE = worker_multiple

##########
# common #
##########

all: $(TARGET_WATER_MULTIPLE)
	 echo "Making water worker for multiple worker setting"

NIMBUS_ROOT = ../..
include ${NIMBUS_ROOT}/Makeinclude

IFLAGS += -I$(LDIR)
H_FILES = $(wildcard *.h)
CPP_FILES = $(wildcard *.cc)
OBJECT_FILES = $(subst .cc,.o,$(CPP_FILES))

PHYSBAM_DIR = $(NIMBUS_ROOT)/application/physbam-lib/Public_Library/
APP_LIB = -lwater_test

####################################
# ENTER application directory here #
####################################

APP_DIR_SINGLE = $(NIMBUS_ROOT)/application/water_test_single
APP_DIR_MULTIPLE = $(NIMBUS_ROOT)/application/water_test_multiple_v1

####################
# ENTER flags here #
####################

LFLAGS_SINGLE = $(LFLAGS)
LDFLAGS_SINGLE = $(LDFLAGS)
IFLAGS_SINGLE = $(IFLAGS)
LFLAGS_SINGLE += -lnimbus -lboost_thread-mt -lboost_system-mt -lpthread $(APP_LIB)
LDFLAGS_SINGLE += -L$(APP_DIR_SINGLE) -Wl,-rpath=$(APP_DIR_SINGLE)
IFLAGS_SINGLE += -I$(APP_DIR_SINGLE) -I$(PHYSBAM_DIR)
USE_SINGLE_FLAG = -D_SINGLE_WORKER

LFLAGS_MULTIPLE = $(LFLAGS)
LDFLAGS_MULTIPLE = $(LDFLAGS)
IFLAGS_MULTIPLE = $(IFLAGS)
LFLAGS_MULTIPLE += -lnimbus -lboost_thread-mt -lboost_system-mt -lpthread $(APP_LIB)
LDFLAGS_MULTIPLE += -L$(APP_DIR_MULTIPLE) -Wl,-rpath=$(APP_DIR_MULTIPLE)
IFLAGS_MULTIPLE += -I$(APP_DIR_MULTIPLE) -I$(PHYSBAM_DIR)

#############################
# add the dependencies here #
#############################

%.o: %.cc
	$(CPP) $(CFLAGS) $(IFLAGS) -c $< -o $@

####################################
# ENTER rules to build target here #
####################################

$(TARGET_WATER_SINGLE): $(CPP_FILES) $(H_FILES) $(OBJECT_FILES)
	$(CPP) $(CFLAGS) $(USE_SINGLE_WORKER) $(IFLAGS_SINGLE) $(LDFLAGS_SINGLE) $(OBJECT_FILES) -o $(TARGET_WATER_SINGLE) $(LFLAGS_SINGLE)

$(TARGET_WATER_MULTIPLE): $(CPP_FILES) $(H_FILES) $(OBJECT_FILES)
	$(CPP) $(CFLAGS) $(IFLAGS_MULTIPLE) $(LDFLAGS_MULTIPLE) $(OBJECT_FILES) -o $(TARGET_WATER_MULTIPLE) $(LFLAGS_MULTIPLE)

clean: clean-files
	\rm -f $(TARGET_WATER_SINGLE) $(TARGET_WATER_MULTIPLE)
