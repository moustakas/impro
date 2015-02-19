OBJ = idl_fivel.o idl_solve.o idl_temden.o idl_pari.o idl_parii.o idl_fivel_wrapper.o
OBJF = im_fivel.f im_input.f im_output.f im_solve.f im_temden.f im_pari.f im_parii.f 

#  Rules for building the object files.

SHL_EXT = so
SHL_LIBS = idl_fivel.$(SHL_EXT)

F_LD_FLAGS = -O -fPIC -shared -ffixed-line-length-132 -fno-automatic -shared

# The following is the default entry point. This section will determine 
# what system we are on, set the correct flags and call this same makefile
# again with the correct flags.

all : 
	@echo "OS type detected: "`uname`
	@case `uname` in \
	   "SunOS") if [  -d /proc ]; then \
			make libs \
			"SHL_EXT=so" \
			"CC=cc" \
			"FC=f77" \
			"C_LD=cc" \
			"F_LD=f77" \
			"CFLAGS=-K pic -DSPARC -O" \
			"C_LD_FLAGS=-O -DSOLARIS" \
			"FFLAGS=-pic -O" \
			"FOBJS=ftn_only_sun.o ftn_strarr_sun.o" \
			"F_LD_FLAGS=-O" \
			"F_LD_POST= -lF77 -lm -lc"; \
	   	    fi  \
		;; \
	   "AIX") make libs \
			"SHL_EXT=a" \
			"CC=cc" \
			"FC=f77" \
			"C_LD=cc" \
			"F_LD=f77" \
			"C_LD_FLAGS=-bM:SRE -bnoentry" \
			"F_LD_FLAGS=-bM:SRE -bnoentry" \
			"EXAMPLE_C2F_EXPORT=-bE:example_c2f.exp" \
			"EXAMPLE1_EXPORT=-bE:example_c2f.exp" \
			"SIMPLEC2F_EXPORT=-bE:simple_c2f.exp" \
			"FORTRAN_EXPORT=-bE:ftn_only.exp" \
			"FOBJS=ftn_only_ibm.o" \
			"F_LD_POST= -lxlf -lxlf90";; \
	   "HP-UX") make libs \
			"SHL_EXT=sl" \
			"CC=cc" \
			"FC=f77" \
			"C_LD=ld" \
			"F_LD=ld" \
			"CFLAGS=+z" \
			"C_LD_FLAGS=-b" \
			"F_LD_FLAGS=-b" \
			"FFLAGS=+z" \
			"FOBJS=ftn_only_hp.o ftn_strarr_hp.o" ;; \
	   "IRIX" ) make libs \
			"SHL_EXT=so" \
			"CC=cc" \
			"FC=@echo \"NOT Implemented for example: f77 \"" \
			"C_LD=ld" \
			"F_LD=@echo \"NOT Implemented for example: f77 \"" \
			"CFLAGS=-KPIC" \
			"C_LD_FLAGS=-shared";; \
	   "IRIX64" ) make libs \
			"SHL_EXT=so" \
			"CC=cc" \
			"FC=@echo \"NOT Implemented for example: f77 \"" \
			"C_LD=ld" \
			"F_LD=@echo \"NOT Implemented for example: f77 \"" \
			"CFLAGS=-KPIC -ansi -32" \
			"C_LD_FLAGS=-shared -32";; \
	   "OSF1" )  make libs \
			"SHL_EXT=so" \
			"CC=cc" \
			"FC=f77" \
			"C_LD=ld" \
			"F_LD=ld" \
			"CFLAGS=-DOSF1" \
			"C_LD_FLAGS=-S -shared" \
			"F_LD_FLAGS=-S -shared" \
			"FOBJS=ftn_only_osf1.o" \
			"F_LD_POST= -lUfor -lfor -lFutil -lm -lots -lc";; \
	   "Linux" ) make libs \
			"SHL_EXT=so" \
			"CC=cc" \
			"F77=g77" \
			"FC=g77" \
			"C_LD=ld" \
			"F_LD=ld" \
			"FFLAGS=-fPIC -O --shared -fno-second-underscore -fugly-assumed -Wuninitialized " \
			"CFLAGS=-fPIC -O -I$(IDL_DIR)/external" \
			"C_LD_FLAGS=-shared -lm";; \
	   *) echo "This system is not supported" ;; \
       esac

im_fivel : 
	g77 -o im_fivel $(OBJF)

libs : $(SHL_LIBS)

idl_fivel.$(SHL_EXT) : $(OBJ)
	$(FC) $(F_LD_FLAGS) -o idl_fivel.$(SHL_EXT) $(OBJ) $(F_LD_POST)

.SUFFIXES: .c .o .f

.c.o :
	$(CC) $(CFLAGS) -c $*.c

.f.o :
	$(F77) $(FFLAGS) -c $*.f

clean: 
	rm *~ *.o *.so
