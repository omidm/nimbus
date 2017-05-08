
TARGET = graph_partitioner

.PHONY: default check clean

default: check $(TARGET)

NIMBUS_ROOT = ../../..
include $(NIMBUS_ROOT)/Makeinclude

CFILES = $(wildcard *.cc)
OBJFILES = $(subst .cc,.o,$(CFILES))

GRAPH_LIB_DIR = $(NIMBUS_ROOT)/applications/graph/graph_library
GRAPH_LIB_CFILES = $(wildcard $(GRAPH_LIB_DIR)/*.cc)
GRAPH_LIB_OBJFILES = $(subst .cc,.o,$(GRAPH_LIB_CFILES))
OBJFILES += $(GRAPH_LIB_OBJFILES)

$(TARGET): $(OBJFILES)
	$(CPP) $(CFLAGS) $(OBJFILES) -o $(TARGET) $(LFLAGS) $(LDFLAGS)

%.o: %.cc
	$(CPP) $(CFLAGS) $(IFLAGS) -c $< -o $@ 

clean: clean-files
	make -C $(GRAPH_LIB_DIR) clean
	\rm -rf $(TARGET)
