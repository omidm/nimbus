
TARGET = bg_process

.PHONY: default check $(TARGET) clean clean-logs

default: check $(TARGET)

NIMBUS_ROOT = ../..
include $(NIMBUS_ROOT)/Makeinclude

CFILES = $(wildcard *.cc)
HFILES = $(wildcard *.h)
OBJFILES = $(subst .cc,.o,$(CFILES))

$(TARGET): $(OBJFILES)
	$(CPP) $(CFLAGS) $(OBJFILES) -o $(TARGET) $(LDFLAGS) $(LFLAGS)

%.o: %.cc
	$(CPP) $(CFLAGS) $(IFLAGS) -c $< -o $@ 

clean: clean-files
	@\rm -f $(TARGET)

clean-logs:
	@\rm -f *log*
	@\rm -f *.txt
	@\rm -rf core

