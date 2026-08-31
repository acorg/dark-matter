import importlib
import sys
from pathlib import Path
from unittest import TestCase

from Bio.Data.CodonTable import TranslationError

sys.dont_write_bytecode = True
sys.path.insert(0, str(Path(__file__).parent.parent / "bin"))
_script = importlib.import_module("analyse-isolate")

classifyChange = _script.classifyChange
classifyMinorVariantBase = _script.classifyMinorVariantBase


class _RaisingTable(dict):
    """Like AmbiguousForwardTable: raises TranslationError for unknown codons."""

    def get(self, key):
        if key not in self:
            raise TranslationError(f"ambiguous codon: {key}")
        return self[key]


def _offsetInfo(featureName="ORF1ab", aaOffset=0):
    return {"featureName": featureName, "reference": {"aaOffset": aaOffset}}


def _clinInfo(aa="A"):
    return {"aa": aa}


def _consInfo(aa="A"):
    return {"aa": aa}


class TestClassifyChange(TestCase):
    """Tests for classifyChange()."""

    def testSynonymousChange(self):
        """Same amino acid on both sides returns 'synonymous' and no furin flag."""
        synInfo, furin = classifyChange(_clinInfo("A"), _consInfo("A"), _offsetInfo())
        self.assertEqual("synonymous", synInfo)
        self.assertEqual("", furin)

    def testNonSynonymousChange(self):
        """Different amino acids produce a non-synonymous description."""
        synInfo, furin = classifyChange(
            _clinInfo("A"),
            _consInfo("V"),
            _offsetInfo(featureName="ORF1ab", aaOffset=100),
        )
        self.assertEqual("non-synonymous (ORF1ab, A101V)", synInfo)
        self.assertEqual("", furin)

    def testNonSynonymousAaOffsetIsOneBased(self):
        """The amino-acid offset displayed is one-based (aaOffset 0 → position 1)."""
        synInfo, _ = classifyChange(
            _clinInfo("M"),
            _consInfo("I"),
            _offsetInfo(featureName="nucleocapsid phosphoprotein", aaOffset=9),
        )
        self.assertEqual("non-synonymous (nucleocapsid phosphoprotein, M10I)", synInfo)

    def testFurinSiteInsideRange(self):
        """Non-synonymous change in the surface glycoprotein furin range is flagged."""
        _, furin = classifyChange(
            _clinInfo("R"),
            _consInfo("A"),
            _offsetInfo(featureName="surface glycoprotein", aaOffset=681),
        )
        self.assertEqual(" FURIN SITE!", furin)

    def testFurinSiteAtLowerBoundary(self):
        """aaOffset == 670 is excluded from the furin range (strict >)."""
        _, furin = classifyChange(
            _clinInfo("R"),
            _consInfo("A"),
            _offsetInfo(featureName="surface glycoprotein", aaOffset=670),
        )
        self.assertEqual("", furin)

    def testFurinSiteAtUpperBoundary(self):
        """aaOffset == 690 is excluded from the furin range (strict <)."""
        _, furin = classifyChange(
            _clinInfo("R"),
            _consInfo("A"),
            _offsetInfo(featureName="surface glycoprotein", aaOffset=690),
        )
        self.assertEqual("", furin)

    def testFurinSiteWrongGene(self):
        """Furin-range offset outside 'surface glycoprotein' gene is not flagged."""
        _, furin = classifyChange(
            _clinInfo("R"),
            _consInfo("A"),
            _offsetInfo(featureName="ORF1ab", aaOffset=681),
        )
        self.assertEqual("", furin)

    def testFurinSiteSynonymousNotFlagged(self):
        """A synonymous change in the furin range is not flagged."""
        synInfo, furin = classifyChange(
            _clinInfo("R"),
            _consInfo("R"),
            _offsetInfo(featureName="surface glycoprotein", aaOffset=681),
        )
        self.assertEqual("synonymous", synInfo)
        self.assertEqual("", furin)


class TestClassifyMinorVariantBase(TestCase):
    """Tests for classifyMinorVariantBase()."""

    _TABLE = {"ATG": "M", "GGG": "G", "AAA": "K", "TTT": "F", "CTG": "L"}

    def testMatchesClinicBase(self):
        """When base equals clinBase the result is 'consensus base'."""
        result = classifyMinorVariantBase(
            base="A",
            clinBase="A",
            clinCod="ATG",
            mvCod="ATG",
            codonTable=self._TABLE,
            offsetInfo=_offsetInfo(),
        )
        self.assertEqual("consensus base", result)

    def testSynonymousVariant(self):
        """Same translated amino acid on both codons → 'synonymous'."""
        result = classifyMinorVariantBase(
            base="G",
            clinBase="T",
            clinCod="ATG",
            mvCod="ATG",
            codonTable=self._TABLE,
            offsetInfo=_offsetInfo(),
        )
        self.assertEqual("synonymous", result)

    def testNonSynonymousVariant(self):
        """Different translated amino acids → formatted AA-change description."""
        result = classifyMinorVariantBase(
            base="G",
            clinBase="T",
            clinCod="ATG",
            mvCod="GGG",
            codonTable=self._TABLE,
            offsetInfo=_offsetInfo(featureName="spike", aaOffset=4),
        )
        self.assertEqual("M5G, spike", result)

    def testNonSynonymousAaOffsetIsOneBased(self):
        """The AA position in the returned string is one-based."""
        result = classifyMinorVariantBase(
            base="G",
            clinBase="T",
            clinCod="ATG",
            mvCod="GGG",
            codonTable=self._TABLE,
            offsetInfo=_offsetInfo(featureName="ORF1ab", aaOffset=0),
        )
        self.assertEqual("M1G, ORF1ab", result)

    def testFeatureNameIncludedInResult(self):
        """The feature name is included verbatim in the non-synonymous description."""
        result = classifyMinorVariantBase(
            base="G",
            clinBase="T",
            clinCod="ATG",
            mvCod="GGG",
            codonTable=self._TABLE,
            offsetInfo=_offsetInfo(
                featureName="nucleocapsid phosphoprotein", aaOffset=0
            ),
        )
        self.assertIn("nucleocapsid phosphoprotein", result)

    def testClinicCodonAmbiguousFallsBackToX(self):
        """When clinCod raises TranslationError, clinCodT is set to 'X'."""
        table = _RaisingTable({"GGG": "G"})
        result = classifyMinorVariantBase(
            base="G",
            clinBase="T",
            clinCod="NNN",
            mvCod="GGG",
            codonTable=table,
            offsetInfo=_offsetInfo(featureName="spike", aaOffset=0),
        )
        self.assertEqual("X1G, spike", result)

    def testMvCodonAmbiguousFallsBackToX(self):
        """When mvCod raises TranslationError, mvCodT is set to 'X'."""
        table = _RaisingTable({"ATG": "M"})
        result = classifyMinorVariantBase(
            base="G",
            clinBase="T",
            clinCod="ATG",
            mvCod="NNN",
            codonTable=table,
            offsetInfo=_offsetInfo(featureName="spike", aaOffset=0),
        )
        self.assertEqual("M1X, spike", result)

    def testBothCodonsAmbiguousGivesSynonymous(self):
        """When both codons raise TranslationError, result is 'synonymous'."""
        table = _RaisingTable()
        result = classifyMinorVariantBase(
            base="G",
            clinBase="T",
            clinCod="NNN",
            mvCod="NAN",
            codonTable=table,
            offsetInfo=_offsetInfo(),
        )
        self.assertEqual("synonymous", result)
