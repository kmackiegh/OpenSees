include ../../../../Makefile.def


OBJS       = ExponentialTS.o ExponentialTS2D.o \
		ElasticTS.o ElasticTS2D.o \
		BilinearTS.o BilinearTS2D.o \
		PDRExponentialTS.o PDRExponentialTS2D.o


all: $(OBJS) 

# Miscellaneous

tidy:
	@$(RM) $(RMFLAGS) Makefile.bak *~ #*# core

clean: tidy
	@$(RM) $(RMFLAGS) $(OBJS) *.o

spotless: clean

wipe: spotless

# DO NOT DELETE THIS LINE -- make depend depends on it.



