
TARGET = graph_partitioner

.PHONY: default check clean

default: check $(TARGET)

NIMBUS_ROOT = ../../..
include $(NIMBUS_ROOT)/Makeinclude

CFILES = $(wildcard *.cc)
OBJFILES = $(subst .cc,.o,$(CFILES))

GRAPH_LIB_DIR = ../graph_library/
GRAPH_LIB_CFILES = $(wildcard $(GRAPH_LIB_DIR)/*.cc)
GRAPH_LIB_OBJFILES = $(subst .cc,.o,$(GRAPH_LIB_CFILES))
OBJFILES += $(GRAPH_LIB_OBJFILES)

$(TARGET): $(OBJFILES)
	$(CPP) $(CFLAGS) $(OBJFILES) -o $(TARGET) $(LFLAGS) $(LDFLAGS)

# graph library objects should not depend on protobuf, so separate rule
$(GRAPH_LIB_DIR)%.o: $(GRAPH_LIB_DIR)%.cc
	$(CPP) $(CFLAGS) $(SHARED_FLAGS) $(IFLAGS) -c $< -o $@

%.o: %.cc
	$(CPP) $(CFLAGS) $(IFLAGS) -c $< -o $@ 

clean: clean-files
	make -C $(GRAPH_LIB_DIR) clean
	\rm -rf $(TARGET)
