import base64
import io
import json
import math
import os
import tempfile

import dash
import dash_bio as dashbio
import numpy as np
import pandas as pd
import plotly.graph_objects as go
from dash import Dash, html, dcc, Output, Input, no_update, State

import modifications
from bam_handler import BAMFileHandler

modifications_dict = {
    "C+m?": "5-Methylcytosine",
    "C+m.": "5-Methylcytosine",
    "C+h?": "5-Hydroxymethylcytosine",
    "C+h.": "5-Hydroxymethylcytosine",
}
reverse_modifications_dict = {v: k for k, v in modifications_dict.items()}

app = Dash(__name__)
app.layout = html.Div([
    html.H1(children='Modiview', style={'textAlign': 'center'}),
    html.Div(id='read-num-span'),
    dcc.Input(
        id='read_number',
        type='number',
        min=1,
        step=1,
        style={'textAlign': 'center', 'width': '5%'},
        value=1),
    html.Div(id='read-num-output-container'),
    "Select modifications:",
    dcc.Dropdown(
        id='mod_dropdown',
        multi=True
    ),
    dashbio.AlignmentChart(
        id='alignment-viewer',
        data='>read1\nGATTACA',
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
        id='upload-file',
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
    ),
    dcc.Store(id='uploaded-file-path'),
    dcc.Store(id='previous-file-path'),
])


@app.callback(
    [Output('read-num-span', 'children'), Output('read_number', 'max')],
    Input('uploaded-file-path', 'data')
)
def update_read_num_span(file_path):
    if file_path is None:
        return 'No file uploaded yet.', no_update
    bam_handler = BAMFileHandler(file_path)
    read_count = bam_handler.get_read_count()
    return f"Enter the read number (1-{read_count}): ", read_count


@app.callback(
    [Output('mod_dropdown', 'options'), Output('mod_dropdown', 'value')],
    [Input('read_number', 'value'), Input('uploaded-file-path', 'data')]
)
def update_mod_dropdown(read_number, file_path):
    if file_path is None:
        return no_update
    bam_handler = BAMFileHandler(file_path)
    read = bam_handler.get_read(read_number)
    mods_short = list(modifications.get_modifications(read).keys())
    mods_full = [modifications_dict[mod] for mod in mods_short]

    return mods_full, mods_full


@app.callback(
    [
        Output('read_number', 'value'),
        Output('zoom-range', 'data')
    ],
    [
        Input('uploaded-file-path', 'data'),
        Input('alignment-viewer', 'eventDatum'),
        Input('read_number', 'value')
    ]
)
def handle_zoom_and_read_number(file_path, alignment_viewer_event, read_number):
    ctx = dash.callback_context

    if not ctx.triggered:
        return no_update
    else:
        input_id = ctx.triggered[0]['prop_id'].split('.')[0]

    if input_id == 'uploaded-file-path':
        return 1, {'start': 0, 'end': 256}
    elif input_id == 'alignment-viewer':
        if alignment_viewer_event is None:
            return no_update
        parsed_event = json.loads(alignment_viewer_event)
        if 'eventType' in parsed_event and parsed_event['eventType'] == 'Zoom':
            start = math.ceil(parsed_event['xStart'])
            end = math.floor(parsed_event['xEnd'])
            bam_handler = BAMFileHandler(file_path)
            read = bam_handler.get_read(read_number)
            sequence = read.query_sequence
            if start < 0 or end > max(len(sequence), 256):
                return no_update
            return no_update, {'start': start, 'end': end}
        else:
            return no_update
    elif input_id == 'read_number':
        return no_update, {'start': 0, 'end': 256}
    else:
        return no_update


@app.callback(
    Output('read-num-output-container', 'children'),
    Input('read_number', 'value'),
    Input('uploaded-file-path', 'data')
)
def update_output(value, file_path):
    if file_path is None:
        return 'File is not uploaded yet.'
    bam_handler = BAMFileHandler(file_path)
    return 'You have selected read number {} with read_id {}'.format(value, bam_handler.get_read_id(value))


@app.callback(
    Output('alignment-viewer', 'data'),
    Input('read_number', 'value'),
    Input('uploaded-file-path', 'data')
)
def update_alignment_chart(value, file_path):
    bam_handler = BAMFileHandler(file_path)
    return bam_handler.get_fasta_sequence(value)


@app.callback(
    Output('modifications-plot', 'figure'),
    Input('zoom-range', 'data'),
    Input('read_number', 'value'),
    Input('mod_dropdown', 'value'),
    Input('uploaded-file-path', 'data')
)
def update_modifications_plot(zoom_range, read_number, mods_selected, file_path):
    if mods_selected is None:
        return no_update
    
    bam_handler = BAMFileHandler(file_path)
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
        width=1500,
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


@app.callback(
    [
        Output('uploaded-file-path', 'data'),
        Output('upload-file', 'contents'),
        Output('upload-file', 'filename')
    ],
    Input('upload-file', 'contents'),
    State('upload-file', 'filename'),
    State('previous-file-path', 'data')
)
def handle_upload(contents, filename, previous_file_path):
    if contents is not None:
        content_type, content_string = contents[0].split(',')
        decoded = base64.b64decode(content_string)
        if 'bam' in filename[0]:
            if previous_file_path is not None and os.path.exists(previous_file_path):
                os.remove(previous_file_path)
            bam_file = io.BytesIO(decoded)
            with tempfile.NamedTemporaryFile(delete=False, suffix=".bam") as fp:
                fp.write(bam_file.read())
                return fp.name, None, None
    return no_update


@app.callback(
    Output('previous-file-path', 'data'),
    Input('uploaded-file-path', 'data')
)
def update_previous_file_path(new_file_path):
    return new_file_path


def generate_modification_list(methylation_positions, nucleotides_shown):
    result = np.zeros(nucleotides_shown)
    for position, probability in methylation_positions:
        result[position] += probability
    return result.tolist()


if __name__ == '__main__':
    app.run(debug=True)
