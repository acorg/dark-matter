def dimensionalIterator(dimensions, maxItems=-1):
    """
    Given a list of n positive integers, return a generator that yields
    n-tuples of coordinates to 'fill' the dimensions. This is like an
    odometer in a car, but the dimensions do not each have to be 10.

    For example: dimensionalIterator((2, 3)) will yield in order
    (0, 0), (0, 1), (0, 2), (1, 0), (1, 1), (1, 2). See the tests in
    test_dimension.py for many more examples.

    A dimension may also be given as '*', to provide a dimension that is
    never exhausted. For example, dimensionalIterator(('*', 2)) yields the
    infinite series (0, 0), (0, 1), (1, 0), (1, 1), (2, 0), (2, 1), ....

    maxItems can be used to limit the number of tuples yielded.
    """
    nDimensions = len(dimensions)
    if nDimensions == 0 or maxItems == 0:
        return
    if any(map(lambda x: x != '*' and x <= 0, dimensions)):
        raise ValueError('Dimensions not all positive! %r' % (dimensions,))
    odometer = [0, ] * nDimensions
    while maxItems != 0:
        yield tuple(odometer)
        maxItems -= 1
        wheel = nDimensions - 1
        while (dimensions[wheel] != '*' and
               odometer[wheel] == dimensions[wheel] - 1 and
               wheel >= 0):
            odometer[wheel] = 0
            wheel -= 1
        if wheel < 0:
            return
        odometer[wheel] += 1
