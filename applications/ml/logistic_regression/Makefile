
TARGET = liblr.so

# add subdirs space separated
SUBDIRS = protobuf_source

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

$(TARGET): $(OBJFILES)
	$(CPP) $(CFLAGS) $(SHARED_FLAGS) $(OBJFILES) -o $(TARGET) $(LDFLAGS) $(LFLAGS)

$(PDIR_COMPILED)%.pb.o: $(PDIR_SOURCE)%.proto
	$(MAKE) -C $(PDIR_SOURCE) ../$@

%.o: %.cc $(PROTO_OBJECT_FILES)
	$(CPP) $(CFLAGS) $(SHARED_FLAGS) $(IFLAGS) -c $< -o $@

clean: clean-files
	for dir in $(SUBDIRS); do \
    $(MAKE) -C $$dir clean; \
  done
	\rm -rf $(TARGET)*

