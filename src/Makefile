
LIB_PATH  = ../lib
LIB_NAME  = libnimbus.so
TARGET    = $(LIB_PATH)/$(LIB_NAME)

default: $(TARGET)

NIMBUS_ROOT = ..
include $(NIMBUS_ROOT)/Makeinclude

# add subdirs space separated
SUBDIRS = shared scheduler worker data application_utils

.PHONY: subdirs $(SUBDIRS) clean

$(SUBDIRS): shared
	$(MAKE) -C $@

SHARED_CFILES = $(wildcard shared/*.cc)
SCHED_CFILES  = $(wildcard scheduler/*.cc)
WORKER_CFILES = $(wildcard worker/*.cc) $(wildcard worker/app_data_managers/*.cc) $(wildcard worker/worker_job_graph/*.cc)
DATA_CFILES   = $(wildcard data/*.cc) $(wildcard data/physbam/*.cc) $(wildcard data/app_data/*.cc)
APP_U_CFILES  = $(wildcard application_utils/*.cc)

CFILES = $(SHARED_CFILES) $(SCHED_CFILES) $(WORKER_CFILES) $(DATA_CFILES) $(APP_U_CFILES) 
OBJFILES = $(subst .cc,.o,$(CFILES))

SHARED_PROTO_OBJECT_FILES = $(wildcard shared/protobuf_compiled/*.pb.o)
OBJFILES += $(SHARED_PROTO_OBJECT_FILES)

DATA_PROTO_OBJECT_FILES = $(wildcard data/physbam/protobuf_compiled/*.pb.o)
OBJFILES += $(DATA_PROTO_OBJECT_FILES)

LFLAGS   = -lboost_thread -lboost_system -lpthread -lprotobuf -lboost_program_options -ldl -lleveldb

# INSTALL_FLAG is not needed for linux
ifdef OS_DARWIN
  INSTALL_FLAG = -install_name @rpath/$(LIB_NAME)
endif

$(TARGET): $(SUBDIRS)
	mkdir -p $(LIB_PATH)
# $(CPP) $(CFLAGS) $(SHARED_FLAGS) $(OBJFILES) -o $(TARGET)
	$(CPP) $(CFLAGS) $(SHARED_FLAGS) $(OBJFILES) -o $(TARGET) $(INSTALL_FLAG) $(LDFLAGS) $(LFLAGS)

clean: clean-files
	for dir in $(SUBDIRS); do \
    $(MAKE) -C $$dir clean; \
  done
	\rm -f $(TARGET)

