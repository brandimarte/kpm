#
#  Makefile for fdf example and fdf module.
#  

# Some systems do not have 'ranlib'. If so,
# use "echo" instead of "ranlib".
RANLIB = ranlib

# Extra flags for library creation by the 'ar' command
# Note that the 'ar' command can itself be specified by
# defining the AR variable. In most 'make' programs, AR is a
# built-in variable.
ARFLAGS_EXTRA =
AR = xiar

default: module

module:  libfdf.a
	cp libfdf.a fdf.mod ..

#
libfdf.a: io_for_fdf_tests.o fdf.o fdf_mod.o parse.o
	$(AR) $(ARFLAGS_EXTRA) cru libfdf.a io_for_fdf_tests.o fdf.o    \
	fdf_mod.o parse.o
	-$(RANLIB) libfdf.a
#
fdf_mod.o: fdf_mod.f parse.o
fdf.o:		fdf.h
#
#---------------------------------------------------
OBJS=	sample.o fdf.o io_for_fdf_tests.o
#
sample: $(OBJS)
	$(FC) $(LDFLAGS) -o sample $(OBJS)
#
check:	sample sample.fdf
	sample < sample.fdf
#
io_sample: io_sample.o io_for_fdf_tests.o
	$(FC) $(LDFLAGS) -o io_sample io_sample.o io.o
#
io_check: io_sample
	io_sample
clean:
	rm -f *.o libfdf.a *.mod fdf_mod.o io_sample.o sample io_sample
	rm -f OUTPUT FDF.debug file2 file3 sample.out
	rm -f loop_prof_*
	rm -f *.pc *.pcl
#
sample.o:	fdfdefs.h



