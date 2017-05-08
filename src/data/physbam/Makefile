
TARGET = physbam

# add subdirs space separated
SUBDIRS = protobuf_source

.PHONY: default check $(TARGET) clean

default: check $(TARGET) $(SUBDIRS)

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

IFLAGS += -I$(NIMBUS_ROOT)/application/physbam-lib/Public_Library

$(TARGET): $(OBJFILES)

$(PDIR_COMPILED)%.pb.o: $(PDIR_SOURCE)%.proto
	$(MAKE) -C $(PDIR_SOURCE) ../$@

%.o: %.cc $(PROTO_OBJECT_FILES)
	$(CPP) $(CFLAGS) $(SHARED_FLAGS) $(IFLAGS) -c $< -o $@

clean: clean-files
	@for dir in $(SUBDIRS); do \
    $(MAKE) -C $$dir clean; \
  done

