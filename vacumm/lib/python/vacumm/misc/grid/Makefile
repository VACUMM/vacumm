PLAT:=$(shell uname -m)
#export LDFLAGS="-Wall -undefined dynamic_lookup -bundle -arch i386 -arch x86_64"
#export FFLAGS="-m32 -m64"

all: $(patsubst %.f90,_%_.so,$(wildcard *.f90))

_%_.so: %.f90
	f2py -c --fcompiler=gnu95  -m $(patsubst %.so,%,$@) $<
