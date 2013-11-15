
all: check application
# all: application

NIMBUS_ROOT = ../..
include $(NIMBUS_ROOT)/Makeinclude

LIBRARY = libstencil_1d_multi.so

CFLAGS += -fPIC
CFILES = $(wildcard *.cc)
HFILES = $(wildcard *.h)
OBJFILES = $(subst .cc,.o,$(CFILES))
LFLAGS += -lboost_thread-mt -lboost_system-mt -lnimbus -lprotobuf
SHARED_FLAGS = -shared -fPIC

PDIR_SOURCE = protobuf_source/
PDIR_COMPILED = protobuf_compiled/
PROTO_FILES = $(wildcard $(PDIR_SOURCE)*.proto)
TEMP_OBJECT_FILES = $(subst .proto,.pb.o,$(PROTO_FILES))
PROTO_OBJECT_FILES = $(subst $(PDIR_SOURCE),$(PDIR_COMPILED),$(TEMP_OBJECT_FILES))

ifdef OS_DARWIN
  LINK_FLAG = -install_name @rpath/$(LIBRARY)
endif

application: $(LIBRARY)

$(LIBRARY): proto_t $(OBJFILES) $(HFILES)
	$(CPP) $(SHARED_FLAGS) $(CFLAGS) $(IFLAGS) $(LDFLAGS) $(LFLAGS) $(PROTO_OBJECT_FILES) $(OBJFILES) -o $(LIBRARY) $(LINK_FLAG)

%.o: %.cc
	$(CPP) $(CFLAGS) $(IFLAGS) -c $< -o $@ 

proto_t:
	cd $(PDIR_SOURCE); make; cd ..

clean: clean-files
	\rm -f $(PDIR_COMPILED)*.pb.*
	\rm -f $(LIBRARY)

