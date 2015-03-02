FC = ifort
FFLAGS =  -align common -r8 -w -i8 -save -g -O2

OBJS=	func.o\
	prepfield.o\
        fit.o\
        ga_carroll.o

fit_ff.exe:	$(OBJS) 
	$(FC) $(FFLAGS) $(OBJS) -o $@



clean:
	rm -f $(OBJS) fit_ff.exe

