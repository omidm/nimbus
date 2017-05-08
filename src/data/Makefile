
TARGET = data

# add subdirs space separated
SUBDIRS = physbam app_data

.PHONY: default check $(TARGET) $(SUBDIRS) clean

default: check $(TARGET) $(SUBDIRS)

NIMBUS_ROOT = ../..
include $(NIMBUS_ROOT)/Makeinclude

CFILES = $(wildcard *.cc)
HFILES = $(wildcard *.h)
OBJFILES = $(subst .cc,.o,$(CFILES))

$(TARGET): $(OBJFILES)

%.o: %.cc
	$(CPP) $(CFLAGS) $(SHARED_FLAGS) $(IFLAGS) -c $< -o $@ 

$(SUBDIRS):
	$(MAKE) -C $@

clean: clean-files
	@for dir in $(SUBDIRS); do \
    $(MAKE) -C $$dir clean; \
  done

