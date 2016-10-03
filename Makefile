EXE     = abpackingangle
OFILES  = abpackingangle.o regression.o matrix.o arrays.o
COPT    = -g -Wall -ansi 
INCLUDE = $(HOME)/include
LIB     = $(HOME)/lib

$(EXE) : $(OFILES)
	$(CC) $(COPT) -L$(LIB) -o $@ $(OFILES) -lbiop -lgen -lm -lxml2

.c.o :
	$(CC) $(COPT) -I$(INCLUDE) -c -o $@ $<

clean :
	\rm -f $(OFILES)

distclean : clean
	\rm -f $(EXE)
