TARGET = worker
APPLICATION = application

default: check $(TARGET) 
all: check $(APPLICATION) $(TARGET) 

NIMBUS_ROOT = ../..
include ${NIMBUS_ROOT}/Makeinclude

IFLAGS += -I$(LDIR)
H_FILES = $(wildcard *.h)
CPP_FILES = $(wildcard *.cc)
OBJECT_FILES = $(subst .cc,.o,$(CPP_FILES))

ADIR = ${NIMBUS_ROOT}/application/page_rank/

LFLAGS += -lnimbus -lboost_thread -lboost_system -lpthread -lpage_rank -lprotobuf -lboost_program_options -lboost_filesystem
LDFLAGS += -L$(ADIR) -Wl,-rpath $(ADIR)

# Add the dependencies here
%.o: %.cc
	$(CPP) $(CFLAGS) $(IFLAGS) -c $< -o $@

$(APPLICATION):
	make -C $(ADIR)

$(TARGET): $(CPP_FILES) $(H_FILES) $(OBJECT_FILES) $(APPLICATION)
	$(CPP) $(CFLAGS) $(IFLAGS) $(LDFLAGS) $(LFLAGS) $(OBJECT_FILES) -o $(TARGET) $(LFLAGS)

clean: clean-files
	\rm -f $(TARGET)

clean-logs:
	\rm -f *log*
	\rm -f *.txt
	\rm -rf _db_*

clean-all: clean-files clean-logs
	\rm -f $(TARGET)
	make -C $(ADIR) clean
