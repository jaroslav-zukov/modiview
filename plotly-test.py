from dash import Dash, html, dcc, callback, Output, Input
import pysam
import dash_bio as dashbio
import plotly.express as px
import pandas as pd
from collections import OrderedDict
import numpy as np

bamfile = pysam.AlignmentFile("calls.bam", "rb", check_sq=False)
read_count = bamfile.count(until_eof=True)
bamfile.close()


def get_read(read_number):
    reads = pysam.AlignmentFile("calls.bam", "rb", check_sq=False)
    counter = 0
    for read in reads:
        counter += 1
        if counter == read_number:
            return read


def get_read_id(read_number):
    return get_read(read_number).query_name


def get_fasta_sequence(read_id):
    result = f">{read_id}\n"
    fastafile = pysam.FastaFile("calls.fasta")
    read = fastafile.fetch(read_id)
    result += read

    return result


def get_modifications(read_number):
    read = get_read(read_number)
    sequence = read.query_sequence
    c_positions = []
    for i in range(len(sequence)):
        if sequence[i] == 'C':
            c_positions.append(i)

    modifications = [e for e in read.get_tag("MM").split(";") if e]
    modification_probabilities_encoded = read.get_tag("ML")
    modification_probabilities = [i / 255 for i in modification_probabilities_encoded]

    modification_relative_position_map = OrderedDict()

    modification_count = 0
    for modification in modifications:
        parts = modification.split(",")
        modification_type = parts[0]
        positions = [int(p) for p in parts[1:]]
        modification_count += len(positions)
        modification_relative_position_map[modification_type] = positions

    modification_position_error_map = {}

    if modification_count == len(modification_probabilities):
        counter = 0
        for modification, relative_position in modification_relative_position_map.items():
            absolute_c_positions = []
            for i in range(len(relative_position)):
                if i == 0:
                    absolute_c_positions.append(relative_position[i])
                else:
                    absolute_c_positions.append(absolute_c_positions[i - 1] + 1 + relative_position[i])

            absolute_positions = [c_positions[i] for i in absolute_c_positions]

            modification_position_error_map[modification] = (
                zip(
                    absolute_positions,
                    modification_probabilities[counter: len(absolute_positions)]
                )
            )

    return modification_position_error_map


app = Dash(__name__)

app.layout = html.Div([
    html.H1(children='Trying visualising stuff', style={'textAlign': 'center'}),
    f"Enter the read number (1-{read_count}): ",
    dcc.Input(
        id='range',
        type='number',
        min=1,
        max=read_count,
        step=1,
        style={'textAlign': 'center', 'width': '5%'},
        value=1),
    html.Div(id='read-num-output-container'),
    dashbio.AlignmentChart(
        id='alignment-viewer',
        data=get_fasta_sequence(get_read_id(1)),
        height=200,
        tilewidth=10,
        tileheight=15,
        textsize=10,
        showconservation=False,
        showconsensus=False,
        showgap=False,
        showid=False,
        colorscale='nucleotide',
        showlabel=False,
    ),
    dcc.Graph(id='modifications-plot'),
    html.Div(id='alignment-chart-width', style={'display': 'none'}),
])


@callback(
    Output('read-num-output-container', 'children'),
    Input('range', 'value'))
def update_output(value):
    return 'You have selected read number {} with read_id {}'.format(value, get_read_id(value))


@callback(
    Output('alignment-viewer', 'data'),
    Input('range', 'value'))
def update_alignment_chart(value):
    read_id = get_read_id(value)
    return get_fasta_sequence(read_id)


@callback(
    Output('alignment-viewer', 'showgap'),
    Input('alignment-viewer', 'eventDatum'))
def update_alignment_chart(eventDatum):
    # print(eventDatum)
    return False


@callback(
    Output('modifications-plot', 'figure'),
    Input('range', 'value'))
def update_modifications_plot(value):
    methylation_positions = get_modifications(value)['C+m?']
    result = np.zeros(256)
    for position, probability in methylation_positions:
        if position < 256:
            result[position] += probability

    x_values = list(range(1, 257))
    y_values = result.tolist()

    df = pd.DataFrame({
        'x': x_values,
        'y': y_values,
        # 'color': y_values  # Use the y-values as the color
    })

    fig = px.bar(df, x='x', y='y',
                 # color='blue',
                 color_continuous_scale=[(0, 'red'), (0.5, 'yellow'), (1, 'green')],
                 height=300)
    fig.update_layout(autosize=False, width=2605, coloraxis_showscale=False, yaxis=dict(range=[0, 1]))

    return fig


if __name__ == '__main__':
    app.run(debug=True)
