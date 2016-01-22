include Makefile.inc

CXXFLAGS = 	$(PRJCXXFLAGS) 
LDLIBS   =  $(PRJLDFLAGS)

VPATH=./test:./src

MAIN = example
MAIN_SRC = example.cc APFELgrid.cc 
MAIN_OBJ = $(MAIN_SRC:%.cc=%.o)

.PHONY: all clean
	
all: $(MAIN)

$(MAIN): $(MAIN_OBJ)

clean:
	-$(RM) -f $(MAIN)
	-$(RM) *.o

force_look:
	true
