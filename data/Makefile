all: check data_t physbam app_data

NIMBUS_ROOT = ..
include $(NIMBUS_ROOT)/Makeinclude

CFLAGS += -fPIC
CFILES = $(wildcard *.cc)
HFILES = $(wildcard *.h)
OBJFILES = $(subst .cc,.o,$(CFILES))

PHYSBAM_CFILES = $(wildcard physbam/*.cc)
PHYSBAM_OFILES = $(subst .cc,.o,$(PHYSBAM_CFILES))

IFLAGS += -I$(NIMBUS_ROOT)/application/physbam-lib/Public_Library

data_t: $(OBJFILES) $(HFILES)

%.o: %.cc
	$(CPP) $(CFLAGS) $(IFLAGS) -c $< -o $@ 

clean: clean-files
	echo "  cleaning physbam/"; cd physbam; make clean; cd ../cache; make clean; cd ..

.PHONY: physbam app_data

physbam: 
	cd physbam && make && cd ..

app_data:
	cd app_data && make && cd ..
