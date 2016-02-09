
all: protobuf

TARGET = protobuf

NIMBUS_ROOT = ../../..
include $(NIMBUS_ROOT)/Makeinclude

CFLAGS += -fPIC

PDIR_SOURCE = ../protobuf_source/
PDIR_COMPILED = ../protobuf_compiled/
PROTO_FILES = $(wildcard $(PDIR_SOURCE)*.proto)
GEN_H = $(subst .proto,.pb.h,$(PROTO_FILES))
TEMP_O = $(subst .proto,.pb.o,$(PROTO_FILES))
PROTO_HEADER_FILES = $(subst $(PDIR_SOURCE),$(PDIR_COMPILED),$(GEN_H))
PROTO_OBJECT_FILES = $(subst $(PDIR_SOURCE),$(PDIR_COMPILED),$(TEMP_O))


$(PDIR_COMPILED)%.pb.h: %.proto
	protoc --cpp_out=$(PDIR_COMPILED) $<

$(PDIR_COMPILED)%.pb.o: $(PROTO_HEADER_FILES)
	$(CPP) $(CFLAGS) $(IFLAGS) -c $(subst .o,.cc,$@) -o $@

$(TARGET): $(PROTO_HEADER_FILES) $(PROTO_OBJECT_FILES)

clean:
	\rm -f $(PDIR_COMPILED)*.pb.*


