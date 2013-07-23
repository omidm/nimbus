TARGET = worker
all: check ${TARGET}

NIMBUS_ROOT = ..
include ../Makeinclude

NIMBUS_APP_ROOT = $(NIMBUS_ROOT)/application/fluid-simulation/
LINKED_APP = $(NIMBUS_APP_ROOT)/libWater.so
LINKED_APP_HEADER = $(wildcard $(NIMBUS_APP_ROOT)/public/*.h)

CFILES = $(wildcard *.cc)
HFILES = $(wildcard *.h)
OBJFILES = $(subst .cc,.o,$(CFILES))
LDFLAGS += -L$(NIMBUS_APP_ROOT)
LFLAGS += -lWater

$(TARGET): $(OBJFILES) $(HFILES) $(LINKED_APP) $(LINKED_APP_HEADER)
	$(CPP) $(CFLAGS) $(IFLAGS) $(LDFLAGS) $(OBJFILES) -o $(TARGET) $(LFLAGS)

%.o: %.cc
	$(CPP) $(CFLAGS) $(IFLAGS) -c $< -o $@ 

clean: clean-files
	\rm -f $(TARGET); rm -f output

run: $(TARGET)
	export LD_LIBRARY_PATH=/home/quhang/Git/cloud/src/nimbus/application/fluid-simulation/; ./worker
	
