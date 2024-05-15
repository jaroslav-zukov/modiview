def calculate_modification_probabilities(read):
    encoded_modification_probabilities = read.get_tag("ML")
    return [i / 255 for i in encoded_modification_probabilities]


def calculate_nucleotide_positions(read, nucleotide):
    return [i for i in range(len(read.query_sequence)) if read.query_sequence[i] == nucleotide]


def get_modification_relative_positions_map(read):
    modifications = [mod for mod in read.get_tag("MM").split(";") if mod]
    modification_relative_positions_map = {}
    for modification in modifications:
        parts = modification.split(",")
        modification_type = parts[0]
        positions = [int(p) for p in parts[1:]]
        modification_relative_positions_map[modification_type] = positions
    return modification_relative_positions_map


def calculate_absolute_modification_positions(relative_positions, nucleotide_positions):
    absolute_nucleotide_positions = []
    for i in range(len(relative_positions)):
        if i == 0:
            absolute_nucleotide_positions.append(relative_positions[i])
        else:
            absolute_nucleotide_positions.append(relative_positions[i - 1] + 1 + relative_positions[i])
    return [nucleotide_positions[i] for i in absolute_nucleotide_positions]


def get_modifications(read):
    modification_relative_position_map = get_modification_relative_positions_map(read)
    modification_probabilities = calculate_modification_probabilities(read)

    modification_position_error_map = {}
    counter = 0
    for modification, relative_position in modification_relative_position_map.items():
        nucleotide_positions = calculate_nucleotide_positions(read, modification[0])
        absolute_positions = calculate_absolute_modification_positions(relative_position, nucleotide_positions)

        modification_position_error_map[modification] = (
            zip(
                absolute_positions,
                modification_probabilities[counter: len(absolute_positions)]
            )
        )

    return modification_position_error_map
