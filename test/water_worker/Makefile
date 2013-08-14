
TARGET_WATER = worker
all: $(TARGET_WATER)

NIMBUS_ROOT = ../..
include ${NIMBUS_ROOT}/Makeinclude

IFLAGS += -I$(LDIR)
H_FILES = $(wildcard *.h)
CPP_FILES = $(wildcard *.cc)
OBJECT_FILES = $(subst .cc,.o,$(CPP_FILES))

APP_DIR = $(NIMBUS_ROOT)/application/water_test
PHYSBAM_DIR = $(NIMBUS_ROOT)/application/physbam-lib/Public_Library/
APP_LIB = -lwater_test
LFLAGS += -lnimbus -lboost_thread-mt -lboost_system-mt -lpthread $(APP_LIB)
LDFLAGS += -L$(APP_DIR) -Wl,-rpath=$(APP_DIR)
IFLAGS += -I$(APP_DIR) -I$(PHYSBAM_DIR)

# Add the dependencies here
%.o: %.cc
	$(CPP) $(CFLAGS) $(IFLAGS) -c $< -o $@

# Link to water application here
$(TARGET_WATER): $(CPP_FILES) $(H_FILES) $(OBJECT_FILES)
	$(CPP) $(CFLAGS) $(IFLAGS) $(LDFLAGS) $(OBJECT_FILES) -o $(TARGET_WATER) $(LFLAGS)

clean: clean-files
	\rm -f $(TARGET)
