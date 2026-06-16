# Ditto for mods (They will be in both lists)
MODS=$(wildcard mod*.f90)
MOD_OBJS= mod_constantes.o mod_pot.o mod_base.o mod_fact.o mod_input.o
TARGET = ifort


FC = gfortran
UNAME_S := $(shell uname -s)
SDKROOT ?= $(shell xcrun --show-sdk-path 2>/dev/null)
ifeq ($(UNAME_S),Darwin)
LIBS = -isysroot $(SDKROOT) -framework Accelerate -lpthread -lm
AS_PREFIX = -B/usr/bin/
else
LIBS = -llapack -liomp5 -lpthread -lm -ldl
AS_PREFIX =
endif
FLFLAGS =-g -O3 -fopenmp -march=native $(AS_PREFIX)

%.o : %.f
	$(FC) -c $(FLFLAGS) $<
%.o: %.f90
	$(FC) -c $(FLFLAGS) $<

MODULES= mod_constantes.o mod_pot.o mod_base.o mod_fact.o mod_input.o
OBJ= yumi.o spline1.o pot.o  tmat.o  imprime2.o  ijohnson.o \
     c0003.o  bes.o  v03.o convergence.o \
     v05.o c0005.o threej0.o threej.o sixj.o ninej.o \
     util.o matinv.o

PROGRAM=yumi.x

default: $(PROGRAM)

$(PROGRAM): | $(MODULES)  $(OBJ)
	$(FC) -o $@ $(MODULES) $(OBJ) $(FLFLAGS) $(LIBS)


debug:
	@echo "OBJS = $(OBJ)"
	@echo "MODS = $(MODULES)"
	@echo "MOD_OBJS = $(MOD_OBJS)"
	@echo "PROGRAM = $(PROGRAM)"
	@echo "PRG_OBJ = $(PRG_OBJ)"

clean:
	rm -rf *.o *.mod *~ $(PROGRAM)

.PHONY: debug default clean
