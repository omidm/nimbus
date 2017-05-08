
# add subdirs space separated
SUBDIRS = page_rank graph_partitioner graph_library

.PHONY: subdirs $(SUBDIRS) clean

subdirs: $(SUBDIRS)

$(SUBDIRS): graph_library
	$(MAKE) -C $@

clean:
	for dir in $(SUBDIRS); do \
    $(MAKE) -C $$dir clean; \
  done

