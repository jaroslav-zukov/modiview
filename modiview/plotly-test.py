import json
import math

import dash
import dash_bio as dashbio
import numpy as np
import pandas as pd
import plotly.express as px
from dash import Dash, html, dcc, callback, Output, Input, no_update

import modifications
from modiview.bam_handler import BAMFileHandler

bam_handler = BAMFileHandler("/Users/jaroslav/Projects/modiview/calls.bam")

read_count = bam_handler.get_read_count()
initial_fasta_sequence = bam_handler.get_fasta_sequence(1)

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
    ),
    dcc.Graph(id='modifications-plot'),
    html.Div(id='alignment-chart-width', style={'display': 'none'}),
])


@callback(
    Output('read-num-output-container', 'children'),
    Input('range', 'value'))
def update_output(value):
    return 'You have selected read number {} with read_id {}'.format(value, bam_handler.get_read_id(value))


@callback(
    Output('alignment-viewer', 'data'),
    Input('range', 'value'))
def update_alignment_chart(value):
    return bam_handler.get_fasta_sequence(value)


@callback(
    Output('modifications-plot', 'figure'),
    Input('alignment-viewer', 'eventDatum'),
    Input('range', 'value'))
def update_modifications_plot(alignment_viewer_event, read_number):
    read = bam_handler.get_read(read_number)
    sequence = read.query_sequence
    trigger_id = dash.callback_context.triggered[0]['prop_id'].split('.')[0]

    start = 0
    end = 256
    nucleotides_shown = max(len(sequence), 256)
    if trigger_id == 'alignment-viewer':
        parsed_event = json.loads(alignment_viewer_event)
        if 'eventType' in parsed_event and parsed_event['eventType'] == 'Zoom':
            start = math.ceil(parsed_event['xStart'])
            end = math.floor(parsed_event['xEnd'])
            # Just not to throw value out of range error
            if start < 0 or end > nucleotides_shown:
                return no_update
        else:
            return no_update

    methylation_positions = modifications.get_modifications(read)['C+m?']
    result = np.zeros(nucleotides_shown)
    for position, probability in methylation_positions:
        result[position] += probability

    x_values = list(range(start + 1, end + 1))
    y_values = result.tolist()[start:end]

    df = pd.DataFrame({
        'x': x_values,
        'y': y_values,
    })

    fig = px.line(df, x='x', y='y', height=300)
    fig.update_layout(autosize=False, width=2605, coloraxis_showscale=False, yaxis=dict(range=[0, 1]))

    return fig


if __name__ == '__main__':
    app.run(debug=True)
