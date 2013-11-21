
#################################
# ENTER application target here #
#################################

TARGET_PROJECTION_MULTIPLE = projection_multiple

##########
# common #
##########

all: $(TARGET_PROJECTION_MULTIPLE)
	 echo "Making projection worker for multiple worker setting"

NIMBUS_ROOT = ../..
include ${NIMBUS_ROOT}/Makeinclude

IFLAGS += -I$(LDIR)
H_FILES = $(wildcard *.h)
CPP_FILES = $(wildcard *.cc)
OBJECT_FILES = $(subst .cc,.o,$(CPP_FILES))

PHYSBAM_DIR = $(NIMBUS_ROOT)/application/physbam-new-lib/Public_Library/
APP_LIB = -lprojection

IFLAGS += -I$(PHYSBAM_DIR)

####################################
# ENTER application directory here #
####################################

APP_DIR_MULTIPLE = $(NIMBUS_ROOT)/application/test_quhang/projection/

####################
# ENTER flags here #
####################

LFLAGS_MULTIPLE = $(LFLAGS)
LDFLAGS_MULTIPLE = $(LDFLAGS)
IFLAGS_MULTIPLE = $(IFLAGS)
LFLAGS_MULTIPLE += -lnimbus -lboost_thread -lboost_system -lpthread -lprotobuf $(APP_LIB)
#LFLAGS_MULTIPLE += -lnimbus -lboost_thread -lboost_system -lpthread $(APP_LIB)
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

$(TARGET_PROJECTION_MULTIPLE): $(CPP_FILES) $(H_FILES) $(OBJECT_FILES)
	$(CPP) $(CFLAGS) $(IFLAGS_MULTIPLE) $(LDFLAGS_MULTIPLE) $(OBJECT_FILES) -o $(TARGET_PROJECTION_MULTIPLE) $(LFLAGS_MULTIPLE)

clean: clean-files
	\rm -f $(TARGET_WATER_SINGLE) $(TARGET_WATER_MULTIPLE)
