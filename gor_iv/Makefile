OBJ = gor4-base.o nrutil.o api.o
HDR = gor4-base.h nrutil.h
CFLAGS = -g -O

all: libgor4.so gor4.so

libgor4.so: $(OBJ)
	gcc -shared -Wl,-install_name,libgor4.so -o libgor4.so -fPIC $(OBJ)

gor4.so: gor4.pyx cgor4.pxd libgor4.so
	CFLAGS=-I. LDFLAGS=-L. python setup.py build_ext -i

nrutil.o: nrutil.h
gor4-base.o: gor4-base.h nrutil.h
lib.o: gor4-base.h nrutil.h

clean:
	rm -fr build gor4.c gor4.so libgor4.so $(OBJ) *~ *.pyc
