#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ(2.59)
AC_INIT([Stacks], [0.9993])
AC_CONFIG_AUX_DIR([config])
AM_INIT_AUTOMAKE([-Wall -Werror foreign])
AM_INIT_AUTOMAKE([subdir-objects])
AC_CONFIG_SRCDIR([src/ustacks.cc])
AC_CONFIG_HEADERS([config.h])

# Checks for programs.
AC_PROG_CXX
AM_PROG_CC_C_O

# Checks for libraries.
AC_CHECK_LIB(gomp, omp_set_num_threads,, echo "Unable to locate OpenMP library")

# Checks for header files.
AC_HEADER_DIRENT
AC_HEADER_STDC
AC_CHECK_HEADERS([float.h limits.h stdlib.h string.h])

# Check for OpenMP parallel execution support
AC_OPENMP

# Checks for typedefs, structures, and compiler characteristics.
AC_HEADER_STDBOOL
AC_C_CONST
AC_TYPE_SIZE_T

# Checks for library functions.
AC_FUNC_CLOSEDIR_VOID
AC_FUNC_MALLOC
AC_FUNC_REALLOC
AC_CHECK_FUNCS([floor memset pow sqrt])

AC_CONFIG_FILES([Makefile])
AC_OUTPUT
