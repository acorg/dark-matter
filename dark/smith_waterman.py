def smith_waterman(seq1, seq2, match=1, mismatch=-1, gap=-1):
    """
    Smith-Waterman algorithm
    Local alignment between two fasta files (nucleotides only)

    @param seq1: C{str} of a sequence that should be aligned
    @param seq2: C{str} of a sequence that should be aligned
    @param match: The C{int} match score.
    @param mismatch: The C{int} mismatch score.
    @param gap: The C{int} gap score.
    @return: A C{list} of 3 strings. The first and third contain
        the input sequences, possibly padded with '-'. The second
        contains '|' where the two sequences match, and ' ' where not.
        If either sequence is of zero length, return a C{list} of three empty
        strings.
    """

    seq1 = seq1.upper()
    seq2 = seq2.upper()

    if mismatch >= 0:
        raise ValueError('Mismatch must be negative')
    if gap >= 0:
        raise ValueError('Gap must be negative')
    # Get sequences
    if len(seq1) == 0 or len(seq2) == 0:
        return ['', '', '']
    for x in seq2:
        if x not in 'ACGT':
            raise ValueError('Invalid DNA nucleotide: "%s"' % x)
    for y in seq1:
        if y not in 'ACGT':
            raise ValueError('Invalid DNA nucleotide: "%s"' % y)
    # Initialisation
    d = {'score': 0, 'pointer': None}
    table = [[d for _ in range(len(seq1)+1)] for _ in range(len(seq2)+1)]
    # Fill
    max_score = 0
    max_x = 0
    max_y = 0

    for x in range(len(seq2)+1)[1:]:
        for y in range(len(seq1)+1)[1:]:
            # Calculate match score
            letter1 = seq1[y-1]
            letter2 = seq2[x-1]
            if letter1 == letter2:
                diagonal_score = table[x-1][y-1]['score'] + match
            else:
                diagonal_score = table[x-1][y-1]['score'] + mismatch
            # Calculate gap scores
            up_score = table[x-1][y]['score'] + gap
            left_score = table[x][y-1]['score'] + gap
            # Choose best score
            if diagonal_score <= 0 and up_score <= 0 and left_score <= 0:
                table[x][y] = {'score': 0, 'pointer': None}
            else:
                if diagonal_score >= up_score:
                    if diagonal_score >= left_score:
                        dia = {'score': diagonal_score, 'pointer': 'diagonal'}
                        table[x][y] = dia
                    else:
                        lef = {'score': left_score, 'pointer': 'left'}
                        table[x][y] = lef
                else:
                    if up_score >= left_score:
                        table[x][y] = {'score': up_score, 'pointer': 'up'}
                    else:
                        lef = {'score': left_score, 'pointer': 'left'}
                        table[x][y] = lef
            # Set max score - is this the best way of getting max score
            # considering how the for loop is iterating through the matrix?
            if table[x][y]['score'] >= max_score:
                max_x = x
                max_y = y
                max_score = table[x][y]['score']

    # Trace-back
    align1 = ''
    align2 = ''
    align = ''
    # printing after local alignment
    current_x = len(seq2)
    current_y = len(seq1)
    if max_x != len(seq2) and max_y != len(seq1):
        if len(seq2)-max_x > len(seq1)-max_y:
            while True:
                align2 += seq2[current_x-1]
                align1 += '-'
                align += ' '
                current_x -= 1
                if current_x - max_x == current_y - max_y:
                    break
        elif len(seq2) - max_x < len(seq1) - max_y:
            while True:
                align1 += seq1[current_y-1]
                align2 += '-'
                align += ' '
                current_y -= 1
                if current_x - max_x == current_y - max_y:
                    break

        while True:
            align1 += seq1[current_y-1]
            align2 += seq2[current_x-1]
            align += ' '
            current_x -= 1
            current_y -= 1
            if max_x == current_x and max_y == current_y:
                break
    elif max_x == len(seq2) and max_y != len(seq1):
        while True:
            align1 += seq1[current_y-1]
            align2 += '-'
            align += ' '
            current_y -= 1
            if current_y == max_y:
                break
    elif max_x != len(seq2) and max_y == len(seq1):
        while True:
            align1 += '-'
            align2 += seq2[current_x-1]
            align += ' '
            current_x -= 1
            if current_x == max_x:
                break
    # traceback
    while True:
        arrow = table[max_x][max_y]['pointer']
        if arrow is None:
            break
        elif arrow == 'diagonal':
            align1 += seq1[max_y-1]
            align2 += seq2[max_x-1]
            letter1a = seq1[max_y-1]
            letter2a = seq2[max_x-1]
            if letter1a == letter2a:
                align += '|'
            else:
                align += ' '
            max_x -= 1
            max_y -= 1
        elif arrow == 'left':
            align1 += seq1[max_y-1]
            align2 += '-'
            align += ' '
            max_y -= 1
        else:
            align1 += '-'
            align2 += seq2[max_x-1]
            align += ' '
            max_x -= 1
    # printing before local alignment
    if max_x != 0 and max_y != 0:
        while True:
            align1 += seq1[max_y-1]
            align2 += seq2[max_x-1]
            align += ' '
            max_x -= 1
            max_y -= 1
            if max_x == 0 or max_y == 0:
                break

        if max_x != 0:
            while True:
                align1 += '-'
                align2 += seq2[max_x-1]
                align += ' '
                max_x -= 1
                if max_x == 0:
                    break
        elif max_y != 0:
            while True:
                align1 += seq1[max_y-1]
                align2 += '-'
                align += ' '
                max_y -= 1
                if max_y == 0:
                    break
    elif max_x == 0 and max_y != 0:
        while True:
            align2 += '-'
            align1 += seq1[max_y-1]
            align += ' '
            max_y -= 1
            if max_y == 0:
                break
    elif max_x != 0 and max_y == 0:
        while True:
            align1 += '-'
            align2 += seq2[max_x - 1]
            align += ' '
            max_x -= 1
            if max_x == 0:
                break

    align1 = align1[::-1]
    align2 = align2[::-1]
    align = align[::-1]

    return [align1, align, align2]


if __name__ == '__main__':
    test = smith_waterman('AGGCTATAGGTTCGTATACCGG',
                          'TATATATAC', match=5, mismatch=-4)
    print '\n'.join(test)
