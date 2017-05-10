# Makefile for StencilProbe
#
#
CC = gcc
COPTFLAGS = $(PAPI) -O3
CLDFLAGS = $(PAPI)

# the line below defines timers.  if not defined, will attempt to automatically
# detect available timers.  See cycle.h.
# should be set to -DHAVE_PAPI or -DHAVEGETTIMEOFDAY or unset.
#TIMER = -DHAVE_PAPI

# the library build is added to link against nimbus applications -omidm
TARGET = libstencilprobe.so
.PHONY: default clean
default: probe $(TARGET)
NIMBUS_ROOT = ../../..
include $(NIMBUS_ROOT)/Makeinclude\
# INSTALL_FLAG is not needed for linux
ifdef OS_DARWIN
  INSTALL_FLAG = -install_name @rpath/$(TARGET)
endif
$(TARGET): util.o probe_heat.o util.h cycle.h common.h
	$(CPP) $(CFLAGS) $(SHARED_FLAGS) util.o probe_heat.o -o $(TARGET) $(INSTALL_FLAG) $(LDFLAGS) $(LFLAGS)
%.o: %.c
	$(CPP) $(CFLAGS) $(SHARED_FLAGS) $(IFLAGS) -c $< -o $@

probe:	main.c util.c run.h probe_heat.c cycle.h
	$(CC) $(COPTFLAGS) $(TIMER) -DRANDOMVALUES main.c util.c probe_heat.c $(CLDFLAGS) -o probe

circqueue_probe:	main.c util.c run.h probe_heat_circqueue.c cycle.h
	$(CC) $(COPTFLAGS) $(TIMER) -DRANDOMVALUES -DCIRCULARQUEUEPROBE main.c util.c probe_heat_circqueue.c $(CLDFLAGS) -o probe

timeskew_probe:	main.c util.c run.h probe_heat_timeskew.c cycle.h
	$(CC) $(COPTFLAGS) $(TIMER) -DRANDOMVALUES main.c util.c probe_heat_timeskew.c $(CLDFLAGS) -o probe

oblivious_probe:	main.c util.c run.h probe_heat_oblivious.c cycle.h
	$(CC) $(COPTFLAGS) $(TIMER) -DRANDOMVALUES main.c util.c probe_heat_oblivious.c $(CLDFLAGS) -o probe

blocked_probe:	main.c util.c probe_heat_blocked.c cycle.h
	$(CC) $(COPTFLAGS) $(TIMER) -DRANDOMVALUES main.c util.c probe_heat_blocked.c $(CLDFLAGS) -o probe

test:	main.c util.c run.h probe_heat.c cycle.h  probe_heat_blocked.c probe_heat_oblivious.c probe_heat_timeskew.c probe_heat_circqueue.c
	$(CC) $(COPTFLAGS) -DSTENCILTEST main.test.c util.c probe_heat.c probe_heat_blocked.c probe_heat_oblivious.c probe_heat_timeskew.c probe_heat_circqueue.c $(CLDFLAGS) -o probe

clean:
	rm -f *.o probe	
	rm -f $(TARGET)
