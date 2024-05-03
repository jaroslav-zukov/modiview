import pysam

# Open the unaligned BAM file
bamfile = pysam.AlignmentFile("calls.bam", "rb", check_sq=False)

read = next(bamfile)
print(read.query_name)
print(read.query_sequence)
print(read.query_qualities)
print(read.get_tag("MM"))
print(read.get_tag("ML"))

query_qualities = read.query_qualities

for i in range(len(query_qualities)):
    print(query_qualities[i])

bamfile.close()

