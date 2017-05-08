
TARGET = libpagerank.so

# add subdirs space separated
SUBDIRS = protobuf_source/ ../graph_library/

.PHONY: default check clean

default: check $(TARGET)

NIMBUS_ROOT = ../../..
include $(NIMBUS_ROOT)/Makeinclude

CFILES = $(wildcard *.cc)
HFILES = $(wildcard *.h)
OBJFILES = $(subst .cc,.o,$(CFILES))

PDIR_SOURCE = protobuf_source/
PDIR_COMPILED = protobuf_compiled/
PROTO_FILES = $(wildcard $(PDIR_SOURCE)*.proto)
TEMP_OBJECT_FILES = $(subst .proto,.pb.o,$(PROTO_FILES))
PROTO_OBJECT_FILES = $(subst $(PDIR_SOURCE),$(PDIR_COMPILED),$(TEMP_OBJECT_FILES))
OBJFILES += $(PROTO_OBJECT_FILES)

GRAPH_LIB_DIR = ../graph_library/
GRAPH_LIB_CFILES = $(wildcard $(GRAPH_LIB_DIR)/*.cc)
GRAPH_LIB_OBJFILES = $(subst .cc,.o,$(GRAPH_LIB_CFILES))
OBJFILES += $(GRAPH_LIB_OBJFILES)

$(TARGET): $(OBJFILES)
	$(CPP) $(CFLAGS) $(SHARED_FLAGS) $(OBJFILES) -o $(TARGET) $(LDFLAGS) $(LFLAGS)

$(PDIR_COMPILED)%.pb.o: $(PDIR_SOURCE)%.proto
	$(MAKE) -C $(PDIR_SOURCE) ../$@

# graph library objects should not depend on protobuf, so separate rule
$(GRAPH_LIB_DIR)%.o: $(GRAPH_LIB_DIR)%.cc
	$(CPP) $(CFLAGS) $(SHARED_FLAGS) $(IFLAGS) -c $< -o $@

%.o: %.cc $(PROTO_OBJECT_FILES)
	$(CPP) $(CFLAGS) $(SHARED_FLAGS) $(IFLAGS) -c $< -o $@

clean: clean-files
	@for dir in $(SUBDIRS); do \
    $(MAKE) -C $$dir clean; \
  done
	@\rm -f $(TARGET)*

