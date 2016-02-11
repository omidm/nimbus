
# add subdirs space separated
SUBDIRS = extern src nodes applications

.PHONY: default $(SUBDIRS) physbam clean clean-hard

default: $(SUBDIRS)

$(SUBDIRS):
	$(MAKE) -C $@

nodes: extern src
applications: extern src

physbam: src
	$(MAKE) -C applications/physbam

# add subdirs space separated
CLEAN_SUBDIRS = src nodes applications
clean:
	for dir in $(CLEAN_SUBDIRS); do \
    $(MAKE) -C $$dir clean; \
  done


# add subdirs space separated
CLEAN_H_SUBDIRS = extern src nodes applications applications/physbam
clean-hard:
	for dir in $(CLEAN_SUBDIRS); do \
    $(MAKE) -C $$dir clean; \
  done

