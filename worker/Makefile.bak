TARGET = worker.so

all: check ${TARGET}

NIMBUS_ROOT = ..
include ../Makeinclude


CFILES = $(wildcard *.cc)
HFILES = $(wildcard *.h)
OBJFILES = $(subst .cc,.o,$(CFILES))
LFLAGS += -lboost_thread-mt -lboost_system-mt
SHARED_FLAGS = -shared -fPIC

$(TARGET): $(OBJFILES) $(HFILES)
	$(CPP) $(CFLAGS) $(SHARED_FLAGS) $(IFLAGS) $(LDFLAGS) $(LFLAGS) $(OBJFILES) -o $(TARGET)

%.o: %.cc
	$(CPP) $(CFLAGS) $(IFLAGS) -c $< -o $@ 

clean: clean-files
	\rm -f $(TARGET)

