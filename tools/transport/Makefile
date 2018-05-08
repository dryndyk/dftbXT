
all: buildwire contacts 

buildwire: buildwire.f90
	gfortran -o buildwire buildwire.f90

contacts: contacts.c
	gcc -o contacts contacts.c

clean:
	rm -f buildwire contacts


