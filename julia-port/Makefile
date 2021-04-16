.SUFFIXES:

mxmlib=mxm.a
FFLAGS=-fdefault-real-8 -fdefault-double-8
OBJDIR=obj
OPT_INCDIR=./
CASEDIR=./

FC=mpif90

L0	=	-O0
L2	=	-O2
L3	=	-O3

FL0	=	$(L0) $(FFLAGS) -g -fbacktrace
FL2	=	$(L2) $(FFLAGS)
FL3	=	$(L3) $(FFLAGS)

#NOBJS = $(patsubst %,$(OBJDIR)/%,$(TMP1))

MXMSRC = mxm_wrapper.f mxm_std.f

objdir: 
	@mkdir $(OBJDIR) 2>/dev/null; cat /dev/null 


mxm :
	$(FC) -c $(FL0) -fPIC -o $(OBJDIR)/mxm_wrapper.o   ./mxm_wrapper.f;
	$(FC) -c $(FL0) -fPIC -o $(OBJDIR)/mxm_std.o   	   ./mxm_std.f;
	$(FC) -shared -fPIC -o mxm.so $(OBJDIR)/mxm_wrapper.o $(OBJDIR)/mxm_std.o -L./blasLapack -lblasLapack

clean :
	rm -rf obj
	mkdir obj
