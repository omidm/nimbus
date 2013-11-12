
TARGET = testr
all: $(TARGET)

NIMBUS_ROOT = ../..
include ${NIMBUS_ROOT}/Makeinclude

IFLAGS += -I$(LDIR)
H_FILES = $(wildcard *.h)
CPP_FILES = $(wildcard *.cc)
OBJECT_FILES = $(subst .cc,.o,$(CPP_FILES))

PDIR = protobufs/
PROTO_FILES = $(wildcard $(PDIR)*.proto)
PROTO_OBJECT_FILES = $(subst .proto,.pb.o,$(PROTO_FILES))

 
LFLAGS += -lnimbus -lboost_thread-mt -lboost_system-mt -lpthread -lprotobuf
LDFLAGS += -L.

# Add the dependencies here
%.o: %.cc
	$(CPP) $(CFLAGS) $(IFLAGS) -c $< -o $@ 

$(TARGET): proto_t $(CPP_FILES) $(H_FILES) $(OBJECT_FILES)
	$(CPP) $(CFLAGS) $(IFLAGS) $(LDFLAGS) $(LFLAGS) $(OBJECT_FILES) $(PROTO_OBJECT_FILES) -o $(TARGET) $(LFLAGS)

proto_t:
	cd $(PDIR); make; cd ..

clean: clean-files
	\rm -f $(TARGET)

