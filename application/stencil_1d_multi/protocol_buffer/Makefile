
all:

TARGET = all

NIMBUS_ROOT = ../../..
include $(NIMBUS_ROOT)/Makeinclude

# iCPP = g++
# CFLAGS = -Wall -g -pg -fPIC
# IFLAGS = -I/usr/local/include # For protobufs

PROTO_FILES = $(wildcard *.proto)
PROTO_OBJECT_FILES = $(subst .proto,.pb.o,$(PROTO_FILES))

%.pb.o: %.proto
	protoc --cpp_out=./ $<
	$(CPP) $(CFLAGS) $(IFLAGS) -c $(subst .o,.cc, $@) -o $@

$(TARGET): $(PROTO_OBJECT_FILES)

clean:
	\rm -f *.pb.*


