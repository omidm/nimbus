
# add subdirs space separated
SUBDIRS = bg_process

.PHONY: subdirs $(SUBDIRS) clean

subdirs: $(SUBDIRS)

$(SUBDIRS):
	$(MAKE) -C $@

clean:
	@for dir in $(SUBDIRS); do \
    $(MAKE) -C $$dir clean; \
  done

