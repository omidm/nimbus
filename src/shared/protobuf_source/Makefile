
TARGET = protobuf

.PHONY: default $(TARGET) clean

default: $(TARGET)

NIMBUS_ROOT = ../../..
include $(NIMBUS_ROOT)/Makeinclude

PDIR_SOURCE = ../protobuf_source/
PDIR_COMPILED = ../protobuf_compiled/
PROTO_FILES = $(wildcard $(PDIR_SOURCE)*.proto)
TEMP = $(subst .proto,.pb.o,$(PROTO_FILES))
PROTO_OBJECT_FILES = $(subst $(PDIR_SOURCE),$(PDIR_COMPILED),$(TEMP))


$(PDIR_COMPILED)%.pb.o: %.proto
	protoc --cpp_out=$(PDIR_COMPILED) $<
	$(CPP) $(CFLAGS) $(SHARED_FLAGS) $(IFLAGS) -c $(subst .o,.cc,$@) -o $@

$(TARGET): $(PROTO_OBJECT_FILES)

clean:
	\rm -f $(PDIR_COMPILED)*.pb.*

