EXE    = cta
OFILES = calculate_torsion_angle.o regression.o matrix.o arrays.o
LIBS   = -lbiop -lgen -lm -lxml2
CC     = gcc -I$(HOME)/include -L$(HOME)/lib -g

$(EXE) : $(OFILES)
	$(CC) -o $@ $(OFILES) $(LIBS)

.c.o :
	$(CC) -c $<

clean :
	\rm -f $(OFILES)

distclean : clean
	\rm -f $(EXE)

