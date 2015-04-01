FC      = ifort
FFLAGS = -cm -w -02
#FFLAGS  = -Wall -Wextra -O3 -fimplicit-none -march=native 
#FFLAGS += -pedantic -fbounds-check -fmax-errors=1 -g
#FFLAGS += $(shell pkg-config --cflags plplotd-f95)
#LDFLAGS = $(shell pkg-config --libs plplotd-f95)
LIBS    =

COMPILE = $(FC) $(FFLAGS)
LINK = $(FC) $(LDFLAGS)

TARGET = createchain       # Name of final executable to produce
OBJS = direction.o sizecalc.o energycalc.o potcalc.o chainmc.o seedrn.o # List of object dependencies

$(TARGET): $(OBJS)
	$(LINK) -o $@ $^ $(LIBS)

%.o:%.f90
	$(COMPILE) -c $<
