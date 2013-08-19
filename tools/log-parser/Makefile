
all: check parser

NIMBUS_ROOT = ../..
include $(NIMBUS_ROOT)/Makeinclude

TARGET = parser

CFLAGS += -fPIC

CFILES = $(wildcard *.cc)
HFILES = $(wildcard *.h)
OBJFILES = $(subst .cc,.o,$(CFILES))

LFLAGS += -lboost_thread-mt -lboost_system-mt
SHARED_FLAGS = -shared -fPIC

ifdef OS_DARWIN
  LINK_FLAG = -install_name @rpath/$(LIBRARY)
endif

parser: $(TARGET)

%.o: %.cc
	$(CPP) $(CFLAGS) $(IFLAGS) -c $< -o $@

$(TARGET): $(CFILES) $(HFILES) $(OBJFILES)
	$(CPP) $(CFLAGS) $(IFLAGS) $(LDFLAGS) $(LFLAGS) $(OBJFILES) -o $(TARGET) $(LINK_FLAG)

clean: clean-files
	\rm -f *.o
	\rm -f $(TARGET)

