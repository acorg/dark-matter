from Bio import SeqIO

def print_len(filename):
	read_number = 0
	total_length = 0
	a_count = 0
	t_count = 0
	c_count = 0
	g_count = 0
	n_count = 0
	for record in SeqIO.parse(open(filename, "rU"), "fasta"):
		total_length += len(record)
		read_number += 1
		for base in record:
			if base == "A":
				a_count += 1
			elif base == "T":
				t_count += 1
			elif base == "C":
				c_count += 1
			elif base == "G":
				g_count += 1
			else:
				n_count += 1
	average_length = total_length / read_number
	average_a = a_count / read_number
	average_t = t_count / read_number
	average_c = c_count / read_number
	average_g = g_count / read_number
	average_n = n_count / read_number
	print "Number of reads: ", read_number
	print "Total length: %s bases" % total_length
	print "The average length: %s bases" % average_length
	print "Total number of A: ", a_count
	print "Total number of T: ", t_count
	print "Total number of C: ", c_count
	print "Total number of G: ", g_count
	print "Total number of N: ", n_count
	print "Average number of A: ", average_a
	print "Average number of T: ", average_t
	print "Average number of C: ", average_c
	print "Average number of G: ", average_g
	print "Average number of N: ", average_n
		

print_len("HCoV-EMC_5_reads.fasta")

