OBJ = gor.o nrutil.o lib.o
HDR = gor.h nrutil.h
CFLAGS = -g

all: libgor.so gor_iv.so

libgor.so: $(OBJ)
	gcc -shared -Wl,-install_name,libgor.so -o libgor.so -fPIC $(OBJ)

gor_iv.so: gor_iv.pyx cgor_iv.pxd libgor.so
	CFLAGS=-I. LDFLAGS=-L. python setup.py build_ext -i

nrutil.o: nrutil.h
gor.o: gor.h nrutil.h
lib.o: gor.h nrutil.h

clean:
	rm -fr build gor_iv.so libgor.so $(OBJ) *~ *.pyc
