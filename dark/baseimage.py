import numpy as np


class BaseImage(object):
    """
    Hold a representation of colored sequence bases in a 2D grid suitable
    for placing onto an alignment graph (see utils.py).
    """

    def __init__(self, xRange, yRange, xScale=1, yScale=1):
        # np.ones is passed a (y range, x range) dimension.
        self.data = np.ones((yRange * yScale + 1, xRange * xScale + 1),
                            dtype=(float, 3))
        self.xScale = xScale
        self.yScale = yScale

    def set(self, x, y, value):
        """
        Set the data at (x, y) to value.
        """
        xBase = int(x) * self.xScale
        yBase = int(y) * self.yScale
        for xOffset in range(self.xScale):
            for yOffset in range(self.yScale):
                self.data[yBase + yOffset, xBase + xOffset] = value
