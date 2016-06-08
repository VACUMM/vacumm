
all: 
	cd ../../../.. && make lib
#all: $(patsubst %.f90,_%_.so,$(wildcard *.f90))

_%_.so: %.f90
	f2py -c --fcompiler=gnu95  -m $(patsubst %.so,%,$@) $<
