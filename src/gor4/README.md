## GOR IV

GOR is an algorithm for predicting secondary structure from an amino acid
sequence. It is described in
[GOR Method for Predicting Protein Secondary Structure from Amino Acid Sequence](http://www.ulb.ac.be/di/map/tlenaert/Home_Tom_Lenaerts/INFO-F-208_files/1996%20Garnier.pdf).

## The original code

The GOR IV code in this directory is based on the C code available from
[http://mig.jouy.inra.fr/?q=en/node/85](http://mig.jouy.inra.fr/?q=en/node/85).

The original C code uses 1-based array offsets. I have not attempted to
change that, or to otherwise clean up the code, tempting as it has been.

## Python interface

I used [cffi](https://cffi.readthedocs.org/en/latest/index.html) to make a
Python interface to GOR IV.

Usage is as follows:

```
from dark.gor4 import GOR4

gor4 = GOR4()
result = gor4.predict('DKATIPSESPFAAAEVADGAIVVDIAKMKYETP')

print('Predicted secondary structure', result['predictions'])
print('Prediction probabilities', result['probabilities'])
```

For detailed usage examples, see the tests in `test/test_gor4.py`.

Note that it is possible to pass alternate known sequence and secondary
structure files to `GOR4.__init__` but this is not done in any of the test
examples.

## Files

The files in this directory are as follows:

* `api.c`: An API library suitable for calling via Cython. This has 3
  functions, `initialize`, `predict` and `finalize`.  Initialize returns a
  `struct` containing pointers to `malloc`d memory that is `free`d in
  `finalize`. To get a secondary structure prediction, `predict` must be
  called.
* `build.py`: Python to tell `cffi` how to build the module, and what functions
  are available.
* `gor4-base.c`: A very slightly modified copy of the original `gor.c` file.  The
  modification is that some `#define` constants have been moved into `gor4-base.h`
  so they can also be used in `api.c`.
* `gor4-base.h`: Some `#define` constants that used to be in `gor4-base.c`, plus
  the definition of the `struct` returned by `initialize`.
* `nrutil.{c,h}`: The original files from the GOR IV distribution.
