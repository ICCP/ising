FC       = gfortran
FFLAGS   = -Wall -Wextra -march=native -O3 -fimplicit-none 
#FFLAGS += -pedantic -fbounds-check -fmax-errors=1 -g
FFLAGS  += $(shell pkg-config --cflags plplot-f95)
LDFLAGS2 = $(shell pkg-config --libs plplot-f95)
LIBS     =

COMPILE = $(FC) $(FFLAGS)
LINK = $(FC) $(LDFLAGS)

TARGET1 = ising.exe          # Name of final executable to produce
OBJS1 = rand_tools.o ising.o # List of object dependencies

TARGET2 = post.exe     # Name of final executable to produce
OBJS2 = plot.o post.o # List of object dependencies

$(TARGET1): $(OBJS1)
	$(LINK) -o $@ $^ $(LIBS)

$(TARGET2): $(OBJS2)
	$(LINK) $(LDFLAGS2) -o $@ $^ $(LIBS)

%.o:%.f95
	$(COMPILE) -c $<

.PHONY:clean
clean:
	rm *.exe *.o *.mod
