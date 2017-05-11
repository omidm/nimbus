
TARGET = libheat.so

LIB_PATH = ../lib/
LIB_NAME = libstencilprobe.so

# add subdirs space separated
SUBDIRS = protobuf_source ${LIB_PATH}

.PHONY: default clean

default: $(TARGET)

NIMBUS_ROOT = ../../../
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

# Link with stencilprobe library. Note that the -L option could be a relative path
# because it is immediately used by the compiler, here. However, the rpath
# option should be absolute path so that when the worker executable loads the
# generated libheat.so, it can find the location of libstencilprobe.so.
LFLAGS  += -l:${LIB_NAME}
LDFLAGS += -L${LIB_PATH}
LDFLAGS += -Wl,-rpath `pwd`/${LIB_PATH}

$(TARGET): $(OBJFILES) ${LIB_PATH}/${LIB_NAME}
	$(CPP) $(CFLAGS) $(SHARED_FLAGS) $(OBJFILES) -o $(TARGET) $(LDFLAGS) $(LFLAGS)

$(PDIR_COMPILED)%.pb.o: $(PDIR_SOURCE)%.proto
	$(MAKE) -C $(PDIR_SOURCE) ../$@

%.o: %.cc $(PROTO_OBJECT_FILES)
	$(CPP) $(CFLAGS) $(SHARED_FLAGS) $(IFLAGS) -c $< -o $@

${LIB_PATH}/${LIB_NAME}:
	$(MAKE) -C ../lib/

clean: clean-files
	for dir in $(SUBDIRS); do \
    $(MAKE) -C $$dir clean; \
  done
	\rm -rf $(TARGET)*

