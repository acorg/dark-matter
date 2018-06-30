# Parameters whose values may vary.

VARIABLE_PARAMS = set([
    'effective_search_space',
    'effective_hsp_length',
    'query',
    'query_length',
    'query_letters'
])


def checkCompatibleParams(initialParams, laterParams):
    """
    Check a later set of BLAST parameters against those originally found.

    @param initialParams: A C{dict} with the originally encountered BLAST
        parameter settings.
    @param laterParams: A C{dict} with BLAST parameter settings encountered
        later.
    @return: A C{str} summary of the parameter differences if the parameter
        sets differ, else C{None}.
    """

    # Note that although the params contains a 'date', its value is empty
    # (as far as I've seen). This could become an issue one day if it
    # becomes non-empty and differs between JSON files that we cat
    # together. In that case we may need to be more specific in our params
    # compatible checking.
    err = []
    for param in initialParams:
        if param in laterParams:
            if (param not in VARIABLE_PARAMS and
                    initialParams[param] != laterParams[param]):
                err.append(
                    '\tParam %r initial value %r differs from '
                    'later value %r' % (param, initialParams[param],
                                        laterParams[param]))
        else:
            err.append('\t%r found in initial parameters, not found '
                       'in later parameters' % param)
    for param in laterParams:
        if param not in initialParams:
            err.append('\t%r found in later parameters, not seen in '
                       'initial parameters' % param)

    return 'Summary of differences:\n%s' % '\n'.join(err) if err else None
