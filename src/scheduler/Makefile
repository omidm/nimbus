all: check scheduler_t

NIMBUS_ROOT = ../..
include $(NIMBUS_ROOT)/Makeinclude

# CFLAGS += -fPIC
CFILES = $(wildcard *.cc)
HFILES = $(wildcard *.h)
OBJFILES = $(subst .cc,.o,$(CFILES))
# LFLAGS += -lboost_thread-mt -lboost_system-mt

scheduler_t: $(OBJFILES) $(HFILES)

%.o: %.cc
	$(CPP) $(CFLAGS) $(SHARED_FLAGS) $(IFLAGS) -c $< -o $@ 

clean: clean-files
