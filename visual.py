import matplotlib.pyplot as plt
from Bio.Seq import Seq
import numpy as np

# DNA sequence
import pysam

# Open the unaligned BAM file
bamfile = pysam.AlignmentFile("calls.bam", "rb", check_sq=False)

read = next(bamfile)

seq = read.query_sequence
query_qualities = read.query_qualities

read_errors = [None] * len(query_qualities)

for i in range(len(query_qualities)):
    # Convert the Phred quality score to an error probability
    read_errors[i] = 10 ** (-query_qualities[i] / 10.0)
    

# mm = read.get_tag("MM")
# ml = read.get_tag("ML")

# Plotting
figure = plt.figure(figsize=(28,4))
axes = plt.axes([0,0,1,1])

span = 1 / len(seq)

# Calculate font size based on figure size and sequence length
dpi = figure.get_dpi()
fig_width = figure.get_figwidth()
font_size = ((fig_width * dpi) / len(seq)) * (72/dpi)

for i in range(len(seq)):
    axes.text(i*span, 0.8, seq[i], 
              fontsize=font_size, 
              color='black', 
              ha='center', 
              va='center',
              bbox=dict(facecolor='none', edgecolor='black', boxstyle='Square, pad=0.15', linewidth=0.5),
              family='monospace'
              )
    
# Create a new Axes object for the bar chart
axes2 = plt.axes([0,0,1,0.2])

# Create the bar chart
axes2.bar(np.arange(len(read_errors))*span, read_errors, color='red', width=span)   

plt.show()

bamfile.close()