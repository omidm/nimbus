
# add subdirs space separated
SUBDIRS = extern src

.PHONY: subdirs $(SUBDIRS) clean clean-src

subdirs: $(SUBDIRS)

$(SUBDIRS):
	$(MAKE) -C $@

clean:
	for dir in $(SUBDIRS); do \
    $(MAKE) -C $$dir clean; \
  done

clean-src:
	$(MAKE) -C src/ clean;

