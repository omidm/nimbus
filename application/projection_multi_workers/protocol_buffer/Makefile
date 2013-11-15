
TARGET = all

CPP = g++
CFLAGS = -Wall -g -pg -fPIC

PROTO_FILES = $(wildcard *.proto)
PROTO_OBJECT_FILES = $(subst .proto,.pb.o,$(PROTO_FILES))

%.pb.o: %.proto
	protoc --cpp_out=./ $<
	$(CPP) $(CFLAGS) -c $(subst .o,.cc, $@) -o $@

$(TARGET): $(PROTO_OBJECT_FILES)

clean:
	\rm -f *.pb.*


