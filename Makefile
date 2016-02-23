
# add subdirs space separated
SUBDIRS = extern src nodes applications

.PHONY: default $(SUBDIRS) physbam clean clean-hard clean-logs

default: $(SUBDIRS)

$(SUBDIRS):
	$(MAKE) -C $@

nodes: extern src
applications: extern src

physbam: src
	$(MAKE) -C applications/physbam


# add subdirs space separated
CLEAN_SUBDIRS = src nodes applications
clean: clean-logs
	for dir in $(CLEAN_SUBDIRS); do \
    $(MAKE) -C $$dir clean; \
  done


# add subdirs space separated
CLEAN_H_SUBDIRS = extern src nodes applications applications/physbam
clean-hard: clean-logs
	for dir in $(CLEAN_H_SUBDIRS); do \
    $(MAKE) -C $$dir clean; \
  done


# add subdirs space separated
CLEAN_L_SUBDIRS = nodes/nimbus_controller nodes/nimbus_worker
clean-logs:
	for dir in $(CLEAN_L_SUBDIRS); do \
    $(MAKE) -C $$dir clean-logs; \
  done
	\rm -rf logs/


.PHONY: test test-stencil test-water

test-stencil: extern src nodes applications
	scripts/test-stencil-basic.sh
	scripts/test-stencil-ft.sh

test-water: extern src nodes applications
	scripts/test-water-basic.sh

test: test-stencil test-water

