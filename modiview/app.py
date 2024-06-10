import json
import math
import os

import dash
import dash_bio as dashbio
import numpy as np
import pandas as pd
import plotly.graph_objects as go
from dash import Dash, html, dcc, callback, Output, Input, no_update

import modifications
from bam_handler import BAMFileHandler

bam_file_relative_path = os.path.join(os.path.dirname(__file__), './../calls.bam')
bam_handler = BAMFileHandler(bam_file_relative_path)

read_count = bam_handler.get_read_count()
initial_fasta_sequence = bam_handler.get_fasta_sequence(1)

modifications_dict = {
    "C+m?": "5-Methylcytosine",
    "C+h?": "5-Hydroxymethylcytosine",
}
reverse_modifications_dict = {v: k for k, v in modifications_dict.items()}

initial_modifications_short = list(modifications.get_modifications(bam_handler.get_read(1)).keys())
initial_modifications_full = [modifications_dict[mod] for mod in initial_modifications_short]

app = Dash(__name__)

app.layout = html.Div([
    html.H1(children='Modiview', style={'textAlign': 'center'}),
    f"Enter the read number (1-{read_count}): ",
    dcc.Input(
        id='read_number',
        type='number',
        min=1,
        max=read_count,
        step=1,
        style={'textAlign': 'center', 'width': '5%'},
        value=1),
    html.Div(id='read-num-output-container'),
    "Select modifications:",
    dcc.Dropdown(
        options=initial_modifications_full,
        value=initial_modifications_full,
        multi=True,
        id='mod_dropdown'
    ),
    dashbio.AlignmentChart(
        id='alignment-viewer',
        data=initial_fasta_sequence,
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
        ticksteps=1000,
    ),
    dcc.Graph(id='modifications-plot'),
    dcc.Store(id='zoom-range', data={'start': 0, 'end': 256}),
    dcc.Upload(
        id='upload-image',
        children=html.Div([
            'Drag and Drop or ',
            html.A('Select Files')
        ]),
        style={
            'width': '99%',
            'height': '60px',
            'lineHeight': '60px',
            'borderWidth': '1px',
            'borderStyle': 'dashed',
            'borderRadius': '5px',
            'textAlign': 'center',
            'margin': '10px'
        },
        multiple=True
    )])


@callback(
    Output('zoom-range', 'data'),
    Input('alignment-viewer', 'eventDatum'),
    Input('read_number', 'value')
)
def update_zoom_range(alignment_viewer_event, read_number):
    if alignment_viewer_event is None:
        return dash.no_update
    parsed_event = json.loads(alignment_viewer_event)
    if 'eventType' in parsed_event and parsed_event['eventType'] == 'Zoom':
        start = math.ceil(parsed_event['xStart'])
        end = math.floor(parsed_event['xEnd'])
        read = bam_handler.get_read(read_number)
        sequence = read.query_sequence
        if start < 0 or end > max(len(sequence), 256):
            return no_update
        return {'start': start, 'end': end}
    else:
        return dash.no_update


@callback(
    Output('read-num-output-container', 'children'),
    Input('read_number', 'value'))
def update_output(value):
    return 'You have selected read number {} with read_id {}'.format(value, bam_handler.get_read_id(value))


@callback(
    Output('alignment-viewer', 'data'),
    Input('read_number', 'value'))
def update_alignment_chart(value):
    return bam_handler.get_fasta_sequence(value)


@callback(
    Output('modifications-plot', 'figure'),
    Input('zoom-range', 'data'),
    Input('read_number', 'value'),
    Input('mod_dropdown', 'value'),
)
def update_modifications_plot(zoom_range, read_number, mods_selected):
    read = bam_handler.get_read(read_number)
    sequence = read.query_sequence

    nucleotides_shown = max(len(sequence), 256)
    start = zoom_range['start']
    end = zoom_range['end']

    fig = go.Figure()
    mods = modifications.get_modifications(read)
    x_values = list(range(start + 1, end + 1))

    for modification in mods_selected:
        y_values = generate_modification_list(
            mods[reverse_modifications_dict[modification]],
            nucleotides_shown
        )[start:end]
        df = pd.DataFrame({
            'x': x_values,
            'y': y_values,
        })
        fig.add_trace(
            go.Scatter(
                x=df['x'],
                y=df['y'],
                mode='lines',
                name=modification
            )
        )

    fig.update_layout(
        autosize=False,
        width=2605,
        height=300,
        yaxis=dict(range=[0, 1]),
        legend=dict(
            orientation="h",
            yanchor="bottom",
            y=1.02,
            xanchor="right",
            x=1
        )
    )

    return fig


def generate_modification_list(methylation_positions, nucleotides_shown):
    result = np.zeros(nucleotides_shown)
    for position, probability in methylation_positions:
        result[position] += probability
    return result.tolist()


if __name__ == '__main__':
    app.run(debug=True)
