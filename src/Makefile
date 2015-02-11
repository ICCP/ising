FC      = gfortran
FFLAGS  = -Wall -Wextra -march=native -O3 -fimplicit-none 
#FFLAGS += -pedantic -fbounds-check -fmax-errors=1 -g
#FFLAGS += $(shell pkg-config --cflags plplotd-f95)
#LDFLAGS = $(shell pkg-config --libs plplotd-f95)
BINDIR := bin

COMPILE = $(FC) $(FFLAGS)
LINK = $(FC) $(LDFLAGS)

TARGET = ising.exe			# Name of final executable to produce
OBJS =  global.o  ising.o	# List of object dependencies

$(TARGET): $(OBJS)
	$(LINK) -o $@ $^ $(LIBS)

%.o:%.f95
	$(COMPILE) -c $<

.PHONY:clean
clean:
	rm *.exe *.o *.mod
