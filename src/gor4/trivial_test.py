#!/usr/bin/env python

from dark.gor4 import GOR4

g = GOR4()
seq = 'DKATIPSESPFAAAEVADGAIVVDIAKMKYETPELHVKVGDTVTWINREA'
print 'sequence   =', seq

result = g.predict('DKATIPSESPFAAAEVADGAIVVDIAKMKYETPELHVKVGDTVTWINREA')
print 'prediction =', result['predictions']
print 'probabilities = ', result['probabilities']
