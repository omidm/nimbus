

TARGET = gradient

.PHONY: default check $(TARGET) run clean

default: check $(TARGET)

NIMBUS_ROOT = ../../../..
include $(NIMBUS_ROOT)/Makeinclude

CFILES = $(wildcard *.cc)
HFILES = $(wildcard *.h)
OBJFILES = $(subst .cc,.o,$(CFILES))

$(TARGET): $(OBJFILES)
	$(CPP) $(CFLAGS) $(IFLAGS) $(LDFLAGS) $(LFLAGS) $(OBJFILES) -o $(TARGET) $(LINK_FLAG) $(LFLAGS)

%.o: %.cc
	$(CPP) $(CFLAGS) $(IFLAGS) -c $< -o $@ 

run: check $(TARGET)
	./$(TARGET) -d 10 -i 20 --sn 0.05

clean: clean-files
	\rm -f $(TARGET)

