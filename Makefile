
# add subdirs space separated
SUBDIRS = extern src nodes applications ec2

.PHONY: default $(SUBDIRS) physbam clean clean-hard clean-logs

default: $(SUBDIRS)

$(SUBDIRS):
	$(MAKE) -C $@

src: extern
ec2: extern src
nodes: extern src
applications: extern src

physbam: src
	$(MAKE) -C applications/physbam


# add subdirs space separated
CLEAN_SUBDIRS = src nodes applications ec2
clean: clean-logs
	for dir in $(CLEAN_SUBDIRS); do \
    $(MAKE) -C $$dir clean; \
  done


# add subdirs space separated
CLEAN_H_SUBDIRS = extern src nodes applications ec2 applications/physbam
clean-hard: clean-logs
	for dir in $(CLEAN_H_SUBDIRS); do \
    $(MAKE) -C $$dir clean; \
  done


# add subdirs space separated
CLEAN_L_SUBDIRS = nodes/nimbus_controller nodes/nimbus_worker ec2/bg_process
clean-logs:
	for dir in $(CLEAN_L_SUBDIRS); do \
    $(MAKE) -C $$dir clean-logs; \
  done
	\rm -rf logs/


.PHONY: stop test test-stencil test-lr test-water


stop:
	scripts/stop-controller.sh
	scripts/stop-workers.sh

test-stencil: extern src nodes applications
	scripts/test-stencil-basic.sh
	scripts/test-stencil-ft.sh

test-lr: extern src nodes applications
	scripts/test-lr-basic.sh
	scripts/test-lr-ft.sh

test-water: extern src nodes applications
	scripts/test-water-basic.sh
	scripts/test-water-ft.sh

test: test-stencil test-lr test-water

