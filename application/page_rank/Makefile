all: check app_t

NIMBUS_ROOT = ../..
include $(NIMBUS_ROOT)/Makeinclude

CFLAGS += -fPIC
LFLAGS += -lnimbus -lboost_thread -lboost_system -lpthread -lprotobuf -lboost_program_options -lboost_filesystem
LDFLAGS += -Wl,-rpath

CFILES = $(wildcard *.cc)
OBJFILES = $(subst .cc,.o,$(CFILES))
TARGET = libpage_rank.so
SHARED_FLAGS = -shared -fPIC

GRAPH_LIB_DIR = ../graph_library
GRAPH_LIB_CFILES = $(wildcard $(GRAPH_LIB_DIR)/*.cc)
GRAPH_LIB_OBJFILES = $(subst .cc,.o,$(GRAPH_LIB_CFILES))

PDIR_SOURCE = protobuf_source/
PDIR_COMPILED = protobuf_compiled/
PROTO_FILES = $(wildcard $(PDIR_SOURCE)*.proto)
TEMP_OBJFILES = $(subst .proto,.pb.o,$(PROTO_FILES))
PROTO_OBJFILES = $(subst $(PDIR_SOURCE),$(PDIR_COMPILED),$(TEMP_OBJFILES))

ifdef OS_DARWIN
  LINK_FLAG = -install_name @rpath/$(LIBRARY)
endif

%.o: %.cc
	$(CPP) $(CFLAGS) $(IFLAGS) -c $< -o $@ 

graph_library_t:
	make -C $(GRAPH_LIB_DIR)

proto_t:
	make -C $(PDIR_SOURCE)

$(TARGET): $(OBJFILES) graph_library_t proto_t
	$(CPP) $(SHARED_FLAGS) $(CFLAGS) $(IFLAGS) $(LDFLAGS) $(LFLAGS) $(PROTO_OBJFILES) $(GRAPH_LIB_OBJFILES) $(OBJFILES) -o $(TARGET) $(LINK_FLAG)

clean: clean-files
	make -C ../graph_library clean
	\rm -f $(PDIR_COMPILED)*.pb.*

app_t: $(TARGET)
