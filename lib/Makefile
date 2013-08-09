
all: check lib

NIMBUS_ROOT = ..
include ../Makeinclude

LIBRARY = libnimbus.so

CFLAGS += -fPIC
CFILES = $(wildcard *.cc)
HFILES = $(wildcard *.h)
OBJFILES = $(subst .cc,.o,$(CFILES))
LFLAGS += -lboost_thread-mt -lboost_system-mt
SHARED_FLAGS = -shared -fPIC

ifdef OS_DARWIN
  LINK_FLAG = -install_name @rpath/$(LIBRARY)
endif

lib: $(LIBRARY)

$(LIBRARY): $(OBJFILES) $(HFILES)
	$(CPP) $(SHARED_FLAGS) $(CFLAGS) $(IFLAGS) $(LDFLAGS) $(LFLAGS) $(OBJFILES) -o $(LIBRARY) $(LINK_FLAG)

%.o: %.cc
	$(CPP) $(CFLAGS) $(IFLAGS) -c $< -o $@ 

clean: clean-files
	\rm -f $(LIBRARY)

