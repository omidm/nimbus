
TARGET = graph_partitioner

.PHONY: default check $(TARGET) clean graph_library_t

default: check $(TARGET)

NIMBUS_ROOT = ../../..
include $(NIMBUS_ROOT)/Makeinclude

CFILES = $(wildcard *.cc)
OBJFILES = $(subst .cc,.o,$(CFILES))

GRAPH_LIB_DIR = $(NIMBUS_ROOT)/applications/graph/graph_library
GRAPH_LIB_CFILES = $(wildcard $(GRAPH_LIB_DIR)/*.cc)
GRAPH_LIB_OBJFILES = $(subst .cc,.o,$(GRAPH_LIB_CFILES))

%.o: %.cc
	$(CPP) $(CFLAGS) $(IFLAGS) -c $< -o $@ 

graph_library_t:
	make -C $(GRAPH_LIB_DIR)

$(TARGET): graph_library_t $(OBJFILES)
	$(CPP) $(CFLAGS) $(GRAPH_LIB_OBJFILES) $(OBJFILES) -o $(TARGET) $(LFLAGS) $(LDFLAGS)

clean: clean-files
	make -C $(GRAPH_LIB_DIR) clean
	\rm -rf $(TARGET)
