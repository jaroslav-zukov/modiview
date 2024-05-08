import pysam

# Open the unaligned BAM file
bamfile = pysam.AlignmentFile("calls.bam", "rb", check_sq=False)

# read = next(bamfile)
# print(read.query_name)
# print(read.query_sequence)
# print(read.query_qualities)
# print(read.get_tag("MM"))
# print(read.get_tag("ML"))

# query_qualities = read.query_qualities

# for i in range(len(query_qualities)):
#     print(query_qualities[i])

print(bamfile.count(until_eof=True))

fastafile = pysam.FastaFile("calls.fasta")

read = fastafile.fetch("6a4800da-842d-4702-9a8e-8bf3bcc264b9")

# Now 'read' contains the sequence of the read with the specified ID
print(read)




fastafile.close()

bamfile.close()

