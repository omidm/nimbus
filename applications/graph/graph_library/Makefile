
TARGET = graph_library_t

.PHONY: default clean $(TARGET)  

default: check $(TARGET)

NIMBUS_ROOT = ../../..
include $(NIMBUS_ROOT)/Makeinclude

CFILES = $(wildcard *.cc)
OBJFILES = $(subst .cc,.o,$(CFILES))

%.o: %.cc
	$(CPP) $(CFLAGS) $(SHARED_FLAGS) $(IFLAGS) -c $< -o $@ 

clean: clean-files

$(TARGET): $(OBJFILES)
