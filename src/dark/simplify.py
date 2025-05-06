def simplifyTitle(title: str, target: str) -> str:
    """
    Simplify a given sequence title.  Given a title, look for the first
    occurrence of target anywhere in any of its words. Return a space-separated
    string of the words of the title up to and including the occurrence of the
    target.  Ignore case.

    E.g.,

      # Suffix
      simplifyTitle('Bovine polyomavirus DNA, complete genome', 'virus') ->
        'Bovine polyomavirus'

      # Prefix
      simplifyTitle('California sea lion polyomavirus 1 CSL6994', 'polyoma') ->
        'California sea lion polyoma'

      # Contained
      simplifyTitle('California sea lion polyomavirus 1 CSL6994', 'yoma') ->
        'California sea lion polyoma'

    title: The string title of the sequence.
    target: The word in the title that we should stop at.

    """
    targetLen = len(target)
    result = []

    for word in title.split():
        if len(word) >= targetLen:
            offset = word.lower().find(target.lower())
            if offset > -1:
                result.append(word[: offset + targetLen])
                break
        result.append(word)
    return " ".join(result)
