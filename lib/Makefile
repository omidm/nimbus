
all: check lib

NIMBUS_ROOT = ..
include ../Makeinclude

TARGET = libnimbus.a

CFILES = $(wildcard *.cc)
HFILES = $(wildcard *.h)
OBJFILES = $(subst .cc,.o,$(CFILES))
LFLAGS += -lboost_thread-mt -lboost_system-mt
SHARED_FLAGS = -shared -fPIC

lib: $(TARGET)

$(TARGET): $(OBJFILES) $(HFILES)
	$(CPP) $(SHARED_FLAGS) $(CPP_FLAGS) $(IFLAGS) $(LDFLAGS) $(LFLAGS) -o $(TARGET) $(OBJFILES) 

%.o: %.cc
	$(CPP) $(CPP_FAGS) $(IFLAGS) -c -o $@ $< $(FLAGS)

clean: clean-files
	\rm -f $(TARGET)

