from dash import Dash, html, dcc, callback, Output, Input
import plotly.express as px
import pandas as pd
import pysam
import dash_bio as dashbio
import urllib.request as urlreq
import json

bamfile = pysam.AlignmentFile("calls.bam", "rb", check_sq=False)
read_count = bamfile.count(until_eof=True)
bamfile.close()

def get_read_id(read_number):
    bamfile = pysam.AlignmentFile("calls.bam", "rb", check_sq=False)
    counter = 0
    for read in bamfile:
        counter += 1
        if counter == read_number:
            return read.query_name
        
data = ">6a4800da-842d-4702-9a8e-8bf3bcc264b9\nGAAAGAAGCTGAGGAGAAGAGCGATTCTGTAGGAAGACCAGCAGTCTCAATTAATATGGACCTCTGAGATCTTTCAAATACTGGACAACCAAACAGACAGCATATACCAGCTGATATGAGGCCCCAAACATACATACAGAAAAGGAATTCCAAGTCTTTGTTCATTCAGGTGCTGTCGATTCCGTTTGTAGTCGTCTGTTTAACCTTAGCAATACGTG"       


app = Dash(__name__)

app.layout = html.Div([
    html.H1(children='Trying visualising stuff', style={'textAlign':'center'}),
    # dcc.Slider(1, read_count, 1, value=1, id='slider-selection', marks=None),
    f"Enter the read number (1-{read_count}): ",
    dcc.Input(id='range', type='number', min=1, max=read_count, step=1, style={'textAlign':'center', 'width': '5%'}, value=1),
    html.Div(id='slider-output-container'),
    dashbio.AlignmentChart(
        id='alignment-viewer',
        data=data,
        height=200,
        tilewidth=10,
        tileheight=15,
        textsize=10,
        showconservation=False,
        showconsensus=False,
        showgap=False,
    ),
    # html.Div(id='alignment-viewer-output')
])

@callback(
    Output('slider-output-container', 'children'),
    Input('range', 'value'))
def update_output(value):
    return 'You have selected read number {} with read_id {}'.format(value, get_read_id(value))


@callback(
    Output('alignment-viewer', 'data'),
    Input('range', 'value'))
def update_alignment_chart(value):
    result = ""
    read_id = get_read_id(value)
    fastafile = pysam.FastaFile("calls.fasta")
    read = fastafile.fetch(read_id)
    result += ">{}\n".format(read_id)
    result += read
    
    return result



if __name__ == '__main__':
    app.run(debug=True)
