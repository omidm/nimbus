
# add subdirs space separated
SUBDIRS = leveldb


.PHONY: subdirs $(SUBDIRS) clean

subdirs: $(SUBDIRS)

all: leveldb_t

$(SUBDIRS):
	$(MAKE) -C $@

clean:
	for dir in $(SUBDIRS); do \
    $(MAKE) -C $$dir clean; \
  done

