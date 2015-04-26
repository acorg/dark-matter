class LocalAlignment(object):
    """
    Smith-Waterman algorithm
    Local alignment between two fasta files (nucleotides only)

    @param seq1: a C{Bio.SeqRecord} instance of a sequence to be aligned
    @param seq2: a C{Bio.SeqRecord} instance of a sequence to be aligned
    @param match: The C{int} match score.
    @param mismatch: The C{int} mismatch score.
    @param gap: The C{int} penalty for opening a gap.
    @param gapExtend: C{int} penalty for extending a gap.
    @param gapExtendDecay: C{float} which decreases the penalty for extending
        a gap.
    @return: A C{list} of strings. For every three lines the first
        and third contain the input sequences, possibly padded with '-'.
        The second contains '|' where the two sequences match,
        and ' ' where not.
        If either sequence is of zero length, raise a ValueError.
        Format of the output is as follows:
        Cigar: (Cigar string)
        Evalue:
        Bitscore:
        Id1 Match start: (int) Match end: (int)
        Id2 Match start: (int) Match end: (int)
        Id1:  1 (seq) 50
        [lines to show matches]
        Id2:  1 (seq) 50
    """

    def __init__(self, seq1, seq2, match=1, mismatch=-1, gap=-1,
                 gapExtend=-1, gapExtendDecay=0.0):
        self.seq1Seq = seq1.seq.upper()
        self.seq1ID = seq1.id
        self.seq2Seq = seq2.seq.upper()
        self.seq2ID = seq2.id
        self.match = match
        self.mismatch = mismatch
        self.gapOpen = gap
        self.gapExtend = gapExtend
        self.gapExtendDecay = gapExtendDecay

        if self.mismatch >= 0:
            raise ValueError('Mismatch must be negative')
        if self.gapOpen >= 0:
            raise ValueError('Gap must be negative')
        if self.gapExtend > 0:
            raise ValueError('Gap extension penalty cannot be positive')

        if len(self.seq1Seq) == 0:
            raise ValueError('Empty sequence: %s' % self.seq1ID)
        if len(self.seq2Seq) == 0:
            raise ValueError('Empty sequence: %s' % self.seq2ID)

        for nt in self.seq1Seq:
            if nt not in 'ACGT':
                raise ValueError('Invalid DNA nucleotide: "%s"' % nt)
        for nt in self.seq2Seq:
            if nt not in 'ACGT':
                raise ValueError('Invalid DNA nucleotide: "%s"' % nt)

    def _initialise(self):
        """
        Initialises table with dictionary.
        """
        d = {'score': 0, 'pointer': None, 'ins': 0, 'del': 0}
        cols = len(self.seq1Seq) + 1
        rows = len(self.seq2Seq) + 1
        table = [[d for _ in range(cols)] for _ in range(rows)]
        return table

    def _fillAndTraceback(self, table):
        """
        Perform Local Alignment according to Smith-Waterman Algorithm.
        Fills the table and then traces back from the highest score.
        NB left = deletion and up = insertion wrt seq1
        """
        # Fill
        max_score = 0
        max_row = 0
        max_col = 0

        for row in range(1, len(self.seq2Seq)+1):
            for col in range(1, len(self.seq1Seq)+1):
                # Calculate match score
                letter1 = self.seq1Seq[col-1]
                letter2 = self.seq2Seq[row-1]
                if letter1 == letter2:
                    diagonal_score = (table[row-1][col-1]['score'] +
                                      self.match)
                else:
                    diagonal_score = (table[row-1][col-1]['score'] +
                                      self.mismatch)

                ins_run = table[row-1][col]['ins']
                del_run = table[row][col-1]['del']

                # Calculate gap scores ensuring extension is not > 0
                if table[row-1][col]['ins'] <= 0:
                    ins_score = table[row-1][col]['score'] + self.gapOpen
                else:
                    if self.gapExtend + ins_run * self.gapExtendDecay <= 0.0:
                        ins_score = (table[row-1][col]['score'] +
                                     self.gapExtend +
                                     ins_run * self.gapExtendDecay)
                    else:
                        ins_score = table[row-1][col]['score']

                if table[row-1][col]['del'] <= 0:
                    del_score = table[row][col-1]['score'] + self.gapOpen
                else:
                    if self.gapExtend + del_run * self.gapExtendDecay <= 0.0:
                        del_score = (table[row][col-1]['score'] +
                                     self.gapExtend +
                                     del_run * self.gapExtendDecay)
                    else:
                        del_score = table[row][col-1]['score']

                # Choose best score
                if diagonal_score <= 0 and ins_score <= 0 and del_score <= 0:
                    table[row][col] = {'score': 0, 'pointer': None, 'ins': 0,
                                       'del': 0}
                else:
                    if diagonal_score >= ins_score:
                        if diagonal_score >= del_score:  # diag lef/up
                            diagonal = {'score': diagonal_score,
                                        'pointer': 'diagonal', 'ins': 0,
                                        'del': 0}
                            table[row][col] = diagonal
                        else:  # lef diag/up
                            deletion = {'score': del_score, 'pointer': 'del',
                                        'ins': 0, 'del': del_run + 1}
                            table[row][col] = deletion
                    else:  # up diag
                        if ins_score >= del_score:  # up diag/lef
                            insertion = {'score': ins_score, 'pointer': 'ins',
                                         'ins': ins_run + 1, 'del': 0}
                            table[row][col] = insertion
                        else:  # lef up diag
                            deletion = {'score': del_score, 'pointer': 'del',
                                        'ins': 0, 'del': del_run + 1}
                            table[row][col] = deletion

                # Set max score - is this the best way of getting max score
                # considering how the for loop iterates through the matrix?
                if table[row][col]['score'] >= max_score:
                    max_row = row
                    max_col = col
                    max_score = table[row][col]['score']

        # Traceback
        indexes = {'max_row': max_row, 'max_col': max_col}
        align1 = ''
        align2 = ''
        align = ''

        current_row = max_row
        current_col = max_col

        while True:
            arrow = table[current_row][current_col]['pointer']
            if arrow is None:
                min_row = current_row + 1
                min_col = current_col + 1
                break
            elif arrow == 'diagonal':
                align1 += self.seq1Seq[current_col-1]
                align2 += self.seq2Seq[current_row-1]
                if self.seq1Seq[current_col-1] == self.seq2Seq[current_row-1]:
                    align += '|'
                else:
                    align += ' '
                current_row -= 1
                current_col -= 1
            elif arrow == 'del':
                align1 += self.seq1Seq[current_col-1]
                align2 += '-'
                align += ' '
                current_col -= 1
            elif arrow == 'ins':
                align1 += '-'
                align2 += self.seq2Seq[current_row-1]
                align += ' '
                current_row -= 1
            else:
                raise ValueError('Invalid pointer: %s' % arrow)

        indexes['min_row'] = min_row
        indexes['min_col'] = min_col
        align1 = align1[::-1]
        align2 = align2[::-1]
        align = align[::-1]

        if len(align1) != len(align2):
            raise ValueError('Lengths of locally aligned sequences'
                             ' different')

        return ([align1, align, align2], indexes)

    def _returnCigarString(self, output):
        """
        Return cigar string of aligned sequences.
        @param output: a C{tup} of strings (align1, align, align2)
        @return: a C{str} containing the cigar string
        Eg with input:
        'GGCCCGCA' and 'GG-CTGCA'
        return 2=1D1=1X3=
        """
        cigar = []
        count = 0
        align1 = output[0]
        align2 = output[2]
        for nt1, nt2 in zip(align1, align2):
            if nt1 == nt2:
                cigar.append('=')
            elif nt1 == '-':
                cigar.append('I')
            elif nt2 == '-':
                cigar.append('D')
            else:
                cigar.append('X')
        # Initially create a list of characters,
        # eg ['=', '=', 'D', '=', 'X', '=', '=', '=']
        cigar.append('*')
        # Append an arbitrary character to ensure parsing below functions
        cigarString = ''
        previousCharacter = ''
        count = 0
        first = True
        for character in cigar:
            if first:
                previousCharacter = character
                count += 1
                first = False
            else:
                if character == previousCharacter:
                    count += 1
                else:
                    cigarString += (str(count) + str(previousCharacter))
                    count = 1
                previousCharacter = character
        return cigarString

    def _formatAlignment(self, output, indexes):
        """
        @param output: a C{tup} of strings (align1, align, align2)
        @param indexes: a C{dict} of positions where the alignment begins/
            ends in the input sequences.
        @return: a C{list} of strings which have been formatted to include
            the ID of the sequences and the positions where the alignment
            begins/ends.
        """
        align1 = ''
        align2 = ''
        align = ''

        if len(self.seq1ID) > len(self.seq2ID):
            diff = len(self.seq1ID) - len(self.seq2ID)
            align1 += self.seq1ID
            align2 += (self.seq2ID + ' ' * diff)
            align += (' ' * len(self.seq1ID))
        elif len(self.seq1ID) < len(self.seq2ID):
            diff = len(self.seq2ID) - len(self.seq1ID)
            align1 += (self.seq1ID + ' ' * diff)
            align2 += self.seq2ID
            align += (' ' * len(self.seq2ID))
        else:
            align1 += (self.seq1ID)
            align2 += (self.seq2ID)
            align += (' ' * len(self.seq1ID))

        if len(str(indexes['min_col'])) > len(str(indexes['min_row'])):
            diff = len(str(indexes['min_col'])) - len(str(indexes['min_row']))
            align1 += (' ' + str(indexes['min_col']) + ' ')
            align2 += (' ' + str(indexes['min_row']) + ' ' * (diff + 1))
            align += (' ' * (2 + len(str(indexes['min_col']))))
        elif len(str(indexes['min_col'])) < len(str(indexes['min_row'])):
            diff = len(str(indexes['min_row'])) - len(str(indexes['min_col']))
            align1 += (' ' + str(indexes['min_col']) + ' ' * (diff + 1))
            align2 += (' ' + str(indexes['min_row']) + ' ')
            align += (' ' * (2 + len(str(indexes['min_row']))))
        else:
            align1 += (' ' + str(indexes['min_col']) + ' ')
            align2 += (' ' + str(indexes['min_row']) + ' ')
            align += (' ' * (2 + len(str(indexes['min_row']))))

        align1 += (output[0] + ' ' + str(indexes['max_col']))
        align += (output[1])
        align2 += (output[2] + ' ' + str(indexes['max_row']))

        return [align1, align, align2]

    def createAlignment(self):
        table = self._initialise()
        alignment = self._fillAndTraceback(table)
        output = alignment[0]
        if output[0] is '' or output[2] is '':
            result = ("\nNo alignment between %s and %s \n"
                      % (self.seq1ID, self.seq2ID))
            return result
        else:
            indexes = alignment[1]
            cigar = self._returnCigarString(output)
            text = self._formatAlignment(output, indexes)
            result = "\n".join(text)
            # bitscore =
            # evalue =
            header = ("\nCigar string of aligned region: %s\n"
                      "%s Match start: %d Match end: %d\n"
                      "%s Match start: %d Match end: %d\n"
                      % (cigar, self.seq1ID, indexes['min_col'],
                         indexes['max_col'], self.seq2ID, indexes['min_row'],
                         indexes['max_row']))
            return header + result
