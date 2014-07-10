def smith_waterman(seq1,seq2):

	"""
	Smith-Waterman algorithm
	Local alignment between two fasta files
	"""

	#Get sequences
	if len(seq1) == 0 or len(seq2) == 0:
		empty = []
		return empty
	#Scoring scheme
	match = 1
	mismatch = -1
	gap = -1
	#Initialisation
	mymatrix = [[{'score' : 0, 'pointer' : 'none'} for x in range(len(seq1)+1)] for x in range(len(seq2)+1)]
	#Fill
	max_score = 0
	max_x = 0
	max_y = 0

	for x in range(len(seq2)+1)[1:]:
		for y in range(len(seq1)+1)[1:]:
			#Calculate match score
			letter1 = seq1[y-1]
			letter2 = seq2[x-1]
			if letter1 == letter2:
				diagonal_score = mymatrix[x-1][y-1]['score'] + match
			else:
				diagonal_score = mymatrix[x-1][y-1]['score'] + mismatch
			#Calculate gap scores
			up_score = mymatrix[x-1][y]['score'] + gap
			left_score = mymatrix[x][y-1]['score'] + gap
			#Choose best score
			if diagonal_score <= 0 and up_score <= 0 and left_score <= 0:
				mymatrix[x][y] = {'score' : 0, 'pointer' : 'none'}
			else:
				if diagonal_score >= up_score:
					if diagonal_score >= left_score:
						mymatrix[x][y] = {'score' : diagonal_score, 'pointer' : 'diagonal'}
					else:
						mymatrix[x][y] = {'score' : left_score, 'pointer' : 'left'}
				else:
					if up_score >= left_score:
						mymatrix[x][y] = {'score' : up_score, 'pointer' : 'up'}
					else:
						mymatrix[x][y] = {'score' : left_score, 'pointer' : 'left'}
			#Set max score - is this the best way of getting max score considering how the for loop is iterating through the matrix?
			if mymatrix[x][y]['score'] >= max_score:
				max_x = x
				max_y = y
				max_score = mymatrix[x][y]['score']

	#print 'max_y is %d, max_x is %d, max_score is %d' % (max_y, max_x, max_score)
	#Trace-back
	align1 = ''
	align2 = ''
	align = ''
	#printing after local alignment
	current_x = len(seq2)
	current_y = len(seq1)
	if (max_x != len(seq2)) and (max_y != len(seq1)):
		if (len(seq2)-max_x) > (len(seq1)-max_y):
			while True:
				align2 += seq2[current_x-1]
				align1 += '-'
				align += ' '
				current_x -= 1
				if (current_x - max_x) == (current_y - max_y):
					break
		elif (len(seq2)-max_x) < (len(seq1)-max_y):
			while True:
				align1 += seq1[current_y-1]
				align2 += '-'
				align += ' '
				current_y -= 1
				if (current_x - max_x) == (current_y - max_y):
					break
		else:
			pass

		while True:
			align1 += seq1[current_y-1]
			align2 += seq2[current_x-1]
			align += ' '
			current_x -= 1
			current_y -= 1
			if max_x == current_x and max_y == current_y:
				break
	elif (max_x == len(seq2)) and (max_y != len(seq1)):
		while True:
			align1 += seq1[current_y-1]
			align2 += '-'
			align += ' '
			current_y -= 1
			if current_y == max_y:
				break
	elif (max_x != len(seq2)) and (max_y == len(seq1)):
		while True:
			align1 += '-'
			align2 += seq2[current_x-1]
			align += ' '
			current_x -= 1
			if current_x == max_x:
				break
	else:
		pass
	#traceback
	while True:
		if mymatrix[max_x][max_y]['pointer'] == 'none':
			break
		elif mymatrix[max_x][max_y]['pointer'] == 'diagonal':
			align1 += seq1[max_y-1]
			align2 += seq2[max_x-1]
			align += '|'
			max_x -= 1
			max_y -= 1
		elif mymatrix[max_x][max_y]['pointer'] == 'left':
			align1 += seq1[max_y-1]
			align2 += '-'
			align += ' '
			max_y -= 1
		else:
			align1 += '-'
			align2 += seq2[max_x-1]
			align += ' '
			max_x -= 1
	#printing before local alignment
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
		else:
			pass
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
	else:
		pass



	#for row in mymatrix:
	#	print row

	align1 = align1[::-1]
	align2 = align2[::-1]
	align = align[::-1]

	result = [align1, align, align2]
	return result

#test = smith_waterman('agtcagtcagtc','cgaatcg')
#print "\n".join(test)