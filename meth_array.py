import numpy as np
from array import array

# The data
data = "C+h?,1,34,1,0,4;C+m?,1,34,1,0,4;"
qualities = array('B', [3, 100, 14, 4, 71, 1, 1, 4, 1, 12])

# Parse the data
segments = data.split(';')
segments = [s for s in segments if s.startswith('C+m')]

# Create a zero array with length 218
result = np.zeros(218)

# For each segment
for segment in segments:
    # Split the segment into parts
    parts = segment.split(',')

    # Get the positions and their probabilities
    positions = [int(p) for p in parts[1:-1]]
    probabilities = qualities[-len(positions):]

    # Convert Phred scores to error probabilities
    error_probabilities = [10 ** (-q / 10) for q in probabilities]

    for i in range(1, len(positions)):
        positions[i] = positions[i-1] + 1 + positions[i]

    # Add the error probabilities to the result array
    for pos, prob in zip(positions, error_probabilities):
        result[pos] += prob


print(result.tolist())