
.PHONY: default boost-t leveldb-t protobuf-t

default: boost-t leveldb-t protobuf-t

NIMBUS_ROOT = ..
include $(NIMBUS_ROOT)/Makeinclude

boost-t: boost
boost:
	./install-boost.sh $(OS_DARWIN)

leveldb-t:
	$(MAKE) -C leveldb

protobuf-t: protobuf
protobuf:
	./install-protobuf.sh


.PHONY: clean-boost clean-leveldb clean-protobuf clean

clean-boost:
	rm -rf boost

clean-leveldb:
	$(MAKE) -C leveldb clean

clean-protobuf:
	rm -rf protobuf

clean: clean-boost clean-leveldb clean-protobuf

