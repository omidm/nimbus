
all: check lib

NIMBUS_ROOT = .
include $(NIMBUS_ROOT)/Makeinclude

LIBRARY = libnimbus.so

CFLAGS += -fPIC

SCHED_CFILES  = $(wildcard scheduler/*.cc)
WORKER_CFILES = $(wildcard worker/*.cc)
DATA_CFILES   = $(wildcard data/*.cc)
SHARED_CFILES = $(wildcard lib/*.cc)
CFILES = $(SCHED_CFILES) $(WORKER_CFILES) $(DATA_CFILES) $(SHARED_CFILES)

HFILES = $(wildcard *.h)
OBJFILES = $(subst .cc,.o,$(CFILES))
LFLAGS += -lboost_thread-mt -lboost_system-mt
SHARED_FLAGS = -shared -fPIC

ifdef OS_DARWIN
  LINK_FLAG = -install_name @rpath/$(LIBRARY)
endif

lib: $(LIBRARY)

.PHONY: scheduler worker data shared
scheduler:  
	cd scheduler; make; cd ..

worker:
	cd worker; make; cd ..

data: 
	cd data; make; cd ..

shared:
	cd lib; make -f Makefile2; cd ..

$(LIBRARY): scheduler worker data shared
	$(CPP) $(SHARED_FLAGS) $(CFLAGS) $(IFLAGS) $(LDFLAGS) $(LFLAGS) $(OBJFILES) -o $(LIBRARY) $(LINK_FLAG)

clean: clean-files
	\rm -f $(LIBRARY)

