
all: extern lib

NIMBUS_ROOT = .
include $(NIMBUS_ROOT)/Makeinclude

.PHONY: extern lib

lib:
	cd src && make -j 12 && cd ..

extern:
	cd extern && make -j 12 && cd ..

clean: clean-files
	\rm -f */*.o */*~ */\#*
	cd src; make clean; cd ..
	cd extern; make clean; cd ..

