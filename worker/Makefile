all: check worker_t worker_job_graph app_data_managers

NIMBUS_ROOT = ..
include $(NIMBUS_ROOT)/Makeinclude

CFLAGS += -fPIC
CFILES = $(wildcard *.cc)
HFILES = $(wildcard *.h)
OBJFILES = $(subst .cc,.o,$(CFILES))

worker_t: $(OBJFILES) $(HFILES)

%.o: %.cc
	$(CPP) $(CFLAGS) $(IFLAGS) -c $< -o $@ 

.PHONY: worker_job_graph app_data_managers

clean: clean-files
	make -C worker_job_graph clean && make -C app_data_managers clean

worker_job_graph:
	make -C worker_job_graph

app_data_managers:
	make -C app_data_managers
