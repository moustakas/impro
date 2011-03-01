      subroutine csp_convol_idl(argc,argv)
      integer*4 argc, argv(*)

c argc = a count of the number of arguements being passed to the routine
c argv = a array of memory pointers

c     j = LOC(argc)           !Obtains number of arguements from argc

c convert IDL parameters to standard FORTRAN pass by reference arguments

      call csp_convol(%VAL(argv(1)),%VAL(argv(2)),%VAL(argv(3)), 
     +     %VAL(argv(4)),%VAL(argv(5)))

      return
      end
