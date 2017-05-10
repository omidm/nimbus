
# add subdirs space separated
SUBDIRS = lib heat

.PHONY: subdirs $(SUBDIRS) clean

subdirs: $(SUBDIRS)

$(SUBDIRS): lib
	$(MAKE) -C $@

clean:
	@for dir in $(SUBDIRS); do \
    $(MAKE) -C $$dir clean; \
  done

