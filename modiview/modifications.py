import numpy as np

modifications_dict = {
    "C+m": "5-Methylcytosine",
    "C+h": "5-Hydroxymethylcytosine",
    "C+21839": "N(4)-methylcytosine",
    "A+a": " 6-Methyladenine",
}
reverse_modifications_dict = {v: k for k, v in modifications_dict.items()}


def generate_modification_list(methylation_positions, nucleotides_shown):
    result = np.zeros(nucleotides_shown)
    for position, probability in methylation_positions:
        result[position] += probability
    return result.tolist()


def list_modifications(read):
    mods_short = list(get_modifications(read).keys())
    return [modifications_dict[mod] for mod in mods_short]


def calculate_modification_probabilities(read):
    encoded_modification_probabilities = read.get_tag("ML")
    return [p / 255 for p in encoded_modification_probabilities]


def calculate_nucleotide_positions(read, nucleotide):
    return [
        i
        for i in range(len(read.query_sequence))
        if read.query_sequence[i] == nucleotide
    ]


def get_modification_relative_positions_map(read):
    modifications = [mod for mod in read.get_tag("MM").split(";") if mod]
    modification_relative_positions_map = {}
    for modification in modifications:
        parts = modification.split(",")
        modification_type = parts[0].rstrip(".?")
        positions = [int(p) for p in parts[1:]]
        modification_relative_positions_map[modification_type] = positions
    return modification_relative_positions_map


def calculate_absolute_modification_positions(relative_positions, nucleotide_positions):
    absolute_nucleotide_positions = []
    for i in range(len(relative_positions)):
        if i == 0:
            absolute_nucleotide_positions.append(relative_positions[i])
        else:
            absolute_nucleotide_positions.append(
                absolute_nucleotide_positions[i - 1] + 1 + relative_positions[i]
            )
    return [nucleotide_positions[i] for i in absolute_nucleotide_positions]


def get_modifications(read):
    modification_relative_positions_map = get_modification_relative_positions_map(read)
    modification_probabilities = calculate_modification_probabilities(read)

    modification_position_error_map = {}
    counter = 0
    for modification, relative_positions in modification_relative_positions_map.items():
        nucleotide_positions = calculate_nucleotide_positions(read, modification[0])
        absolute_positions = calculate_absolute_modification_positions(
            relative_positions, nucleotide_positions
        )

        modification_position_error_map[modification] = zip(
            absolute_positions,
            modification_probabilities[counter : counter + len(absolute_positions)],
        )
        counter += len(absolute_positions)

    return modification_position_error_map


def generate_plot_data(modification, read, nucleotides_shown):
    mods = get_modifications(read)
    return generate_modification_list(
        mods[reverse_modifications_dict[modification]], nucleotides_shown
    )
