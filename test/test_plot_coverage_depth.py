import importlib.util
from pathlib import Path
from unittest import TestCase

# Load the hyphen-named script as a module.
_script = Path(__file__).parent.parent / "bin" / "plot-coverage-depth.py"
_spec = importlib.util.spec_from_file_location("plot_coverage_depth", _script)
_mod = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_mod)

assignRegionRows = _mod.assignRegionRows
parseRegion = _mod.parseRegion


class TestParseRegion(TestCase):
    """
    Tests for parseRegion.
    """

    def testStartEnd(self):
        """
        Only start and end given; label and color must be None.
        """
        self.assertEqual((10, 20, None, None), parseRegion("10:20"))

    def testStartEndLabel(self):
        """
        Start, end, and label given; color must be None.
        """
        self.assertEqual((10, 20, "MyGene", None), parseRegion("10:20:MyGene"))

    def testStartEndLabelColor(self):
        """
        All four fields given.
        """
        self.assertEqual((10, 20, "MyGene", "red"), parseRegion("10:20:MyGene:red"))

    def testStartMustBePositive(self):
        """
        A start value of 0 must raise ValueError.
        """
        self.assertRaises(ValueError, parseRegion, "0:10")

    def testEndMustBeGeStart(self):
        """
        An end value less than start must raise ValueError.
        """
        self.assertRaises(ValueError, parseRegion, "10:5")

    def testNonIntegerStart(self):
        """
        A non-integer start must raise ValueError.
        """
        self.assertRaises(ValueError, parseRegion, "abc:10")

    def testMissingEnd(self):
        """
        A region string with only one field must raise ValueError.
        """
        self.assertRaises(ValueError, parseRegion, "10")


class TestAssignRegionRows(TestCase):
    """
    Tests for assignRegionRows.
    """

    def _region(self, start, end):
        return (start, end, None, None)

    def testEmptyRegions(self):
        """
        No regions must return an empty list.
        """
        self.assertEqual([], assignRegionRows([]))

    def testSingleRegion(self):
        """
        A single region must be assigned to row 0.
        """
        self.assertEqual([0], assignRegionRows([self._region(1, 100)]))

    def testNonOverlappingRegions(self):
        """
        Non-overlapping regions must all fit in row 0.
        """
        regions = [self._region(1, 50), self._region(51, 100), self._region(101, 200)]
        self.assertEqual([0, 0, 0], assignRegionRows(regions))

    def testTwoOverlappingRegions(self):
        """
        Two overlapping regions must be placed in rows 0 and 1.
        """
        regions = [self._region(1, 100), self._region(50, 150)]
        self.assertEqual([0, 1], assignRegionRows(regions))

    def testThreeOverlappingRegions(self):
        """
        Three mutually overlapping regions must each get their own row.
        """
        regions = [self._region(1, 100), self._region(10, 90), self._region(20, 80)]
        self.assertEqual([0, 1, 2], assignRegionRows(regions))

    def testAdjacentRegionReusesRow(self):
        """
        A region starting strictly after the previous one ends must reuse row 0.
        """
        regions = [self._region(1, 50), self._region(51, 100)]
        self.assertEqual([0, 0], assignRegionRows(regions))

    def testOutOfOrderNonOverlappingRegions(self):
        """
        Regions given in reverse x-order that do not overlap must all be
        placed in row 0 once sorted by start position.
        """
        # Supplied largest-first; without sorting this would spill into extra rows.
        regions = [self._region(200, 300), self._region(100, 150), self._region(1, 50)]
        self.assertEqual([0, 0, 0], assignRegionRows(regions))

    def testOutOfOrderWithOverlap(self):
        """
        Out-of-order regions where some overlap must still stack correctly.
        """
        # L1=2950-2960, L2=2965-2975, L3=3150-3160 (none overlap)
        # C1=1900-1920, C2=2350-2360 (none overlap with each other or L*)
        regions = [
            self._region(2950, 2960),  # L1
            self._region(2965, 2975),  # L2
            self._region(3150, 3160),  # L3
            self._region(1900, 1920),  # C1
            self._region(2350, 2360),  # C2
        ]
        # None of these overlap, so all should be in row 0.
        self.assertEqual([0, 0, 0, 0, 0], assignRegionRows(regions))
