OBJ = idl_fivel.o idl_solve.o idl_temden.o idl_pari.o idl_parii.o idl_fivel_wrapper.o
OBJF = im_fivel.f im_input.f im_output.f im_solve.f im_temden.f im_pari.f im_parii.f 

#  Rules for building the object files.

SHL_EXT = so
SHL_LIBS = idl_fivel.$(SHL_EXT)

F_LD_FLAGS = -O -fPIC -shared -ffixed-line-length-132 -fno-automatic -shared

im_fivel : 
	gfortran -o im_fivel $(OBJF)

libs : $(SHL_LIBS)

idl_fivel.$(SHL_EXT) : $(OBJ)
	$(FC) $(F_LD_FLAGS) -o idl_fivel.$(SHL_EXT) $(OBJ) $(F_LD_POST)

.SUFFIXES: .c .o .f

.c.o :
	$(CC) $(CFLAGS) -c $*.c

.f.o :
	$(F77) $(FFLAGS) -c $*.f

clean: 
	rm -f *~ *.o *.so im_fivel
