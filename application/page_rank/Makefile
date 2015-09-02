all: check app_t

NIMBUS_ROOT = ../..
include $(NIMBUS_ROOT)/Makeinclude

CFLAGS += -fPIC
LFLAGS += -lnimbus -lboost_thread -lboost_system -lpthread -lprotobuf -lboost_program_options -lboost_filesystem
LDFLAGS += -Wl,-rpath

CFILES = $(wildcard *.cc)
OBJFILES = $(subst .cc,.o,$(CFILES))
TARGET = libpage_rank.so
SHARED_FLAGS = -shared -fPIC

ifdef OS_DARWIN
  LINK_FLAG = -install_name @rpath/$(LIBRARY)
endif

%.o: %.cc
	$(CPP) $(CFLAGS) $(IFLAGS) -c $< -o $@ 

graph_library_t:
	make -C ../graph_library

$(TARGET): $(OBJFILES) graph_library_t
	$(CPP) $(SHARED_FLAGS) $(CFLAGS) $(IFLAGS) $(LDFLAGS) $(LFLAGS) $(OBJFILES) -o $(TARGET) $(LINK_FLAG)

clean: clean-files
	make -C ../graph_library clean

app_t: $(TARGET)
