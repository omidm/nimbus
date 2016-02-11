
TARGET = libnimbus.so

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

$(TARGET): $(SUBDIRS)
	$(CPP) $(CFLAGS) $(SHARED_FLAGS) $(OBJFILES) -o $(TARGET)

clean: clean-files
	for dir in $(SUBDIRS); do \
    $(MAKE) -C $$dir clean; \
  done
	\rm -f $(TARGET)

