
TARGET = app_data_managers

.PHONY: default check $(TARGET) clean

default: check $(TARGET)

NIMBUS_ROOT = ../../..
include $(NIMBUS_ROOT)/Makeinclude

CFILES = $(wildcard *.cc)
HFILES = $(wildcard *.h)
OBJFILES = $(subst .cc,.o,$(CFILES))

$(TARGET): $(OBJFILES)

%.o: %.cc
	$(CPP) $(CFLAGS) $(SHARED_FLAGS) $(IFLAGS) -c $< -o $@ 

clean: clean-files

