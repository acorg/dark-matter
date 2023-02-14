from unittest import TestCase

from dark.graphics import plotAAProperties
from dark.reads import AARead


class TestPlotAAProperties(TestCase):
    """
    Tests for the plotAAProperties function in graphics.py
    """

    # Note that plotAAProperties relies on dark.aa.propertiesForSequence
    # which is tested thoroughly in test_aa.py

    def testNoProperties(self):
        """
        plotAAProperties must run correctly when called with no properties.
        """
        read = AARead("id", "AI")
        self.assertEqual({}, plotAAProperties(read, [], showFigure=False))

    def testOneProperty(self):
        """
        plotAAProperties must run correctly when called with one property.
        """
        read = AARead("id", "AI")
        self.assertEqual(
            {
                "composition": [-1.0, -1.0],
            },
            plotAAProperties(read, ["composition"], showFigure=False),
        )

    def testTwoProperties(self):
        """
        plotAAProperties must run correctly when called with two properties.
        """
        read = AARead("id", "AI")
        self.assertEqual(
            {
                "composition": [-1.0, -1.0],
                "hydropathy": [0.4, 1.0],
            },
            plotAAProperties(read, ["composition", "hydropathy"], showFigure=False),
        )
