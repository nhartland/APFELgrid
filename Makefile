include Makefile.inc

CXXFLAGS = 	$(PRJCXXFLAGS) 
LDLIBS   =  $(PRJLDFLAGS)

VPATH=./test:./src

MAIN = test0
MAIN_SRC = test0.cc 
MAIN_OBJ = $(MAIN_SRC:%.cc=%.o)

.PHONY: all clean
	
all: $(MAIN)

$(MAIN): $(MAIN_OBJ)

clean:
	-$(RM) -f $(MAIN)
	-$(RM) *.o

force_look:
	true
