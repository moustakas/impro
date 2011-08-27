      subroutine idl_fivel_wrapper(argc,argv)

      integer*4 argc, argv(*)

c argc - a count of the number of arguements being passed to the routine
c argv - an array of memory pointers

c convert IDL parameters to standard FORTRAN pass-by-reference arguments

      call idl_fivel(%val(argv(1)), %val(argv(2)), %val(argv(3)), 
     &     %val(argv(4)), %val(argv(5)), %val(argv(6)), %val(argv(7)),
     &     %val(argv(8)), %val(argv(9)), %val(argv(10)), %val(argv(11)), 
     &     %val(argv(12)), %val(argv(13)), %val(argv(14)), 
     &     %val(argv(15)))

      return
      end
