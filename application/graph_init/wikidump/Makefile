all: check graph_initialize_t

NIMBUS_ROOT = ../../..
include $(NIMBUS_ROOT)/Makeinclude

CFLAGS += -fPIC
LFLAGS += -lnimbus -lboost_thread -lboost_system -lpthread -lprotobuf -lboost_program_options -lboost_filesystem
LDFLAGS += -Wl,-rpath

CFILES = $(wildcard *.cc)
OBJFILES = $(subst .cc,.o,$(CFILES))
TARGET = graph_initialize

GRAPH_LIB_DIR = $(NIMBUS_ROOT)/application/graph_library
GRAPH_LIB_CFILES = $(wildcard $(GRAPH_LIB_DIR)/*.cc)
GRAPH_LIB_OBJFILES = $(subst .cc,.o,$(GRAPH_LIB_CFILES))

%.o: %.cc
	$(CPP) $(CFLAGS) $(IFLAGS) -c $< -o $@ 

graph_library_t:
	make -C $(GRAPH_LIB_DIR)

graph_initialize_t: graph_library_t $(OBJFILES)
	$(CPP) $(CFLAGS) $(IFLAGS) $(LDFLAGS) $(LFLAGS) $(GRAPH_LIB_OBJFILES) $(OBJFILES) -o $(TARGET) $(LFLAGS)

clean: clean-files
	\rm -rf $(TARGET)
