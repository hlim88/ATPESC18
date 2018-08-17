AC_PREREQ(2.63)
AC_INIT([MLIFE],
        1.0,
        [wgropp@illinois.edu],
        [mlife],
        [http://www.cs.illinois.edu/~wgropp/advmpi])

AC_CONFIG_HEADER(code-examples-c/mlifeconf.h)

AC_ARG_WITH(mpi,[--with-mpi=dir - location of MPI compilation and run scripts],,with_mpi=/usr/local/bin)
AC_ARG_WITH(pnetcdf,[--with-pnetcdf=dir - location of pNetCDF libraries],,with_pnetcdf=no)
AC_ARG_WITH(fpmpi2,[--with-fpmpi2=dir - location of FPMPI2 libraries],,with_fpmpi2=no)

# Set the C and Fortran FLAGS so that autoconf won't set them to
# its defaults, which aren't appropriate for these examples
: ${CFLAGS=""}
: ${FCFLAGS=""}
AC_SUBST(FCFLAGS)
#
# Try to find the compiler
if test -d "$with_mpi" -a -x "$with_mpi/mpicc" ; then
    MPI_DIR=$with_mpi
fi

# Find the compilers, if they aren't set with the environment variables
if test -z "$MPICC" ; then
    AC_MSG_CHECKING([for how to compile MPI program in C])
    if test -x $MPI_DIR/mpicc ; then
        MPICC=$MPI_DIR/mpicc
    elif test -x $MPI_DIR/bin/mpicc ; then
        MPICC=$MPI_DIR/bin/mpicc
    elif test -x $MPI_DIR/mpcc ; then
        MPICC=$MPIDIR/mpcc
    elif test -x $MPI_DIR/bin/mpcc ; then
        MPICC=$MPIDIR/bin/mpcc
    else
        saveIFS=$IFS
        IFS=":"
	for dir in $PATH ; do
	    if test -x $dir/mpicc ; then
                MPICC=$dir/mpicc
		break
            elif test -x $dir/mpcc ; then
                MPICC=$dir/mpcc
		break
            fi
	done
	IFS=$saveIFS
    fi
    if test -n "$MPICC" ; then
        AC_MSG_RESULT($MPICC)
    else
        AC_MSG_RESULT([cannot determine!])
    fi
fi

if test -z "$MPIFORT" ; then
    AC_MSG_CHECKING([for how to compile MPI program in Fortran])
    if test -x $MPI_DIR/mpifort ; then
        MPIFORT=$MPI_DIR/mpifort
    elif test -x $MPI_DIR/bin/mpif90 ; then
        MPIFORT=$MPI_DIR/bin/mpiff90
    elif test -x $MPI_DIR/mpxlf ; then
        MPIFORT=$MPIDIR/mpxlf
    elif test -x $MPI_DIR/bin/mpxlf90 ; then
        MPIFORT=$MPIDIR/bin/mpxlf90
    else
        saveIFS=$IFS
        IFS=":"
	for dir in $PATH ; do
	    if test -x $dir/mpifort ; then
                MPIFORT=$dir/mpifort
		break
            elif test -x $dir/mpif90 ; then
                MPIFORT=$dir/mpif90
		break
            elif test -x $dir/mpxlf ; then
                MPIFORT=$dir/mpxlf
		break
            fi
	done
	IFS=$saveIFS
    fi
    if test -n "$MPIFORT" ; then
        AC_MSG_RESULT($MPIFORT)
    else
        AC_MSG_RESULT([cannot determine!])
    fi
fi

AC_SUBST(MPICC)
AC_SUBST(MPIFORT)

# Now, make sure that we have something that we can use
if test -n "$MPICC" ; then
    # Note that this may make CC a cross compiler
    CC=$MPICC
else
    AC_MSG_ERROR([You must specify the directory containing the MPI compilation scripts using the --with-mpi=dir option])
fi

AC_PROG_CC
AC_C_VOLATILE
AC_C_RESTRICT
AC_C_INLINE
#
# Find the MPI location
# Try to find drand48
# Try to find nanosleep
AC_CHECK_FUNCS(drand48 nanosleep)
# If we didn't find them, try with a few options
if test "$ac_cv_func_drand48" != yes ; then
    saveCFLAGS=$CFLAGS
    CFLAGS="$CFLAGS -D_XOPEN_SOURCE"
    AC_MSG_CHECKING([for drand48 defined when _XOPEN_SOURCE is defined])
    AC_TRY_LINK([#include <stdlib.h>],[double a=drand48();],pac_cv_func_drand48_with_xopen=yes,pac_cv_func_drand48_with_xopen=no)
    AC_MSG_RESULT($pac_cv_drand48_with_xopen)
    if test "$pac_cv_drand48_with_xopen" != yes ; then
        CFLAGS=$saveCFLAGS
        AC_MSG_WARN([Some examples require the random number generator drand48])
    fi
fi

if test "$with_pnetcdf" != no -a -d "$with_pnetcdf" ; then
    PNETCDF_DIR=$with_pnetcdf
fi

if test "$with_fpmpi2" != no -a -d "$with_fpmpi2" ; then
    FPMPI2_DIR=$with_fpmpi2
fi

AC_SUBST(MPI_DIR)
AC_SUBST(PNETCDF_DIR)
AC_SUBST(FPMPI2_DIR)
AC_OUTPUT([Makefile
	   code-examples-c/Makefile
	   code-examples-f9x/Makefile
	])
