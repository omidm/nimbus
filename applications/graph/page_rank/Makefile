
TARGET = libpagerank.so

PDIR_SOURCE = protobuf_source/
GRAPH_LIB_DIR = ../graph_library

# add subdirs space separated
SUBDIRS = $(PDIR_SOURCE) $(GRAPH_LIB_DIR)

.PHONY: default check $(TARGET) $(SUBDIRS) clean

default: check $(TARGET) $(SUBDIRS)

NIMBUS_ROOT = ../../..
include $(NIMBUS_ROOT)/Makeinclude

CFILES = $(wildcard *.cc)
HFILES = $(wildcard *.h)
OBJFILES = $(subst .cc,.o,$(CFILES))

PDIR_COMPILED = protobuf_compiled/
PROTO_FILES = $(wildcard $(PDIR_SOURCE)*.proto)
TEMP_OBJECT_FILES = $(subst .proto,.pb.o,$(PROTO_FILES))
PROTO_OBJECT_FILES = $(subst $(PDIR_SOURCE),$(PDIR_COMPILED),$(TEMP_OBJECT_FILES))
OBJFILES += $(PROTO_OBJECT_FILES)

GRAPH_LIB_CFILES = $(wildcard $(GRAPH_LIB_DIR)/*.cc)
GRAPH_LIB_OBJFILES = $(subst .cc,.o,$(GRAPH_LIB_CFILES))
OBJFILES += $(GRAPH_LIB_OBJFILES)

$(TARGET): $(SUBDIRS) $(OBJFILES)
	$(CPP) $(CFLAGS) $(SHARED_FLAGS) $(OBJFILES) -o $(TARGET) $(LDFLAGS) $(LFLAGS)

%.o: %.cc
	$(CPP) $(CFLAGS) $(SHARED_FLAGS) $(IFLAGS) -c $< -o $@ 

%.o: $(SUBDIRS)

$(SUBDIRS):
	$(MAKE) -C $@

clean: clean-files
	for dir in $(SUBDIRS); do \
    $(MAKE) -C $$dir clean; \
  done
	\rm -f $(TARGET)

