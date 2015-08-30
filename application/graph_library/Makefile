all: check graph_library_t graph_initialize_t

NIMBUS_ROOT = ../..
include $(NIMBUS_ROOT)/Makeinclude

CFLAGS += -fPIC
LFLAGS += -lnimbus -lboost_thread -lboost_system -lpthread -lprotobuf -lboost_program_options -lboost_filesystem
LDFLAGS += -Wl,-rpath

CFILES = $(wildcard *.cc)
OBJFILES = $(subst .cc,.o,$(CFILES))
TARGET = graph_initialize

%.o: %.cc
	$(CPP) $(CFLAGS) $(IFLAGS) -c $< -o $@ 

$(TARGET): graph_library_t
	$(CPP) $(CFLAGS) $(IFLAGS) $(LDFLAGS) $(LFLAGS) $(OBJFILES) -o $(TARGET) $(LFLAGS)

clean: clean-files

graph_initialize_t: $(TARGET)
graph_library_t: $(OBJFILES)
