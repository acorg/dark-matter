import string
from os.path import basename


def numericallySortFilenames(names):
    """
    Sort (ascending) a list of file names by their numerical prefixes.
    The number sorted on is the numeric prefix of the basename of
    the given filename. E.g., '../output/1.json.bz2' will sort before
    '../output/10.json.bz2'.

    @param: A C{list} of file names, each of whose basename starts with a
        string of digits.
    @return: The sorted C{list} of full file names.
    """

    def numericPrefix(name):
        """
        Find any numeric prefix of C{name} and return it as an C{int}.

        @param: A C{str} file name, whose name possibly starts with digits.
        @return: The C{int} number at the start of the name, else 0 if
            there are no leading digits in the name.
        """
        count = 0
        for ch in name:
            if ch in string.digits:
                count += 1
            else:
                break
        return 0 if count == 0 else int(name[0:count])

    return sorted(names, key=lambda name: numericPrefix(basename(name)))
