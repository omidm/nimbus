
TARGET = worker
all: check $(TARGET)

NIMBUS_ROOT = ../..
include ${NIMBUS_ROOT}/Makeinclude

IFLAGS += -I$(LDIR)
H_FILES = $(wildcard *.h)
CPP_FILES = $(wildcard *.cc)
OBJECT_FILES = $(subst .cc,.o,$(CPP_FILES))

ADIR = ../../application/1d-stencil/
APP_H_FILES = $(wildcard $(ADIR)*.h)
APP_CPP_FILES = $(wildcard $(ADIR)*.cc)
APP_OBJECT_FILES = $(subst .cc,.o,$(APP_CPP_FILES))
 
LFLAGS += -lnimbus -lboost_thread-mt -lboost_system-mt -lpthread -l1d_stencil
LDFLAGS += -L$(ADIR) -Wl,-rpath=$(ADIR)

# Add the dependencies here
%.o: %.cc
	$(CPP) $(CFLAGS) $(IFLAGS) -c $< -o $@

$(ADIR)%.o: $(ADIR)%.cc
	$(CPP) $(CFLAGS) $(IFLAGS) -c $< -o $@

$(TARGET): $(CPP_FILES) $(H_FILES) $(OBJECT_FILES) $(APP_OBJECT_FILES)
	$(CPP) $(CFLAGS) $(IFLAGS) $(LDFLAGS) $(LFLAGS) $(OBJECT_FILES) $(APP_OBJECT_FILES) -o $(TARGET) $(LFLAGS)

clean: clean-files
	\rm -f $(TARGET)

