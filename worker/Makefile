all: check worker_t worker_job_graph

NIMBUS_ROOT = ..
include $(NIMBUS_ROOT)/Makeinclude

CFLAGS += -fPIC
CFILES = $(wildcard *.cc)
HFILES = $(wildcard *.h)
OBJFILES = $(subst .cc,.o,$(CFILES))

worker_t: $(OBJFILES) $(HFILES)

%.o: %.cc
	$(CPP) $(CFLAGS) $(IFLAGS) -c $< -o $@ 

.PHONY: worker_job_graph

clean: clean-files
	cd worker_job_graph; make clean; cd ..

worker_job_graph:
	cd worker_job_graph; make; cd ..
