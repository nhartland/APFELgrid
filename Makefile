include Makefile.inc

CXXFLAGS = 	$(PRJCXXFLAGS) 
LDLIBS   =  $(PRJLDFLAGS)

VPATH=./test:./src

GEN = example_gen
GEN_SRC = example_gen.cc APFELgrid.cc fastkernel.cc fkgenerator.cc
GEN_OBJ = $(GEN_SRC:%.cc=%.o)

CON = example_conv
CON_SRC = example_conv.cc APFELgrid.cc fastkernel.cc fkgenerator.cc
CON_OBJ = $(CON_SRC:%.cc=%.o)

.PHONY: all clean
	
all: $(GEN) $(CON)

$(GEN): $(GEN_OBJ)
$(CON): $(CON_OBJ)

clean:
	-$(RM) -f $(GEN) $(CON)
	-$(RM) *.o

force_look:
	true
