class DatabaseDuplicationError(Exception):
    """An attempt was made to insert a duplicate into a database."""


class NoSuchGenomeError(Exception):
    """A genome accession number could not be found."""


class NoSuchProteinError(Exception):
    """A protein accession number could not be found."""


class ReadLengthsNotIdenticalError(Exception):
    """Not all read lengths are identical."""
