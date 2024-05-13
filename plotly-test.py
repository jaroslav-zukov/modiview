from dash import Dash, html, dcc, callback, Output, Input
import pysam
import dash_bio as dashbio
import plotly.express as px
import pandas as pd

bamfile = pysam.AlignmentFile("calls.bam", "rb", check_sq=False)
read_count = bamfile.count(until_eof=True)
bamfile.close()


def get_read_id(read_number):
    reads = pysam.AlignmentFile("calls.bam", "rb", check_sq=False)
    counter = 0
    for read in reads:
        counter += 1
        if counter == read_number:
            return read.query_name


def get_fasta_sequence(read_id):
    result = f">{read_id}\n"
    fastafile = pysam.FastaFile("calls.fasta")
    read = fastafile.fetch(read_id)
    result += read

    return result


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
    x_values = list(range(1, 257))
    y_values = list(range(1, 257))

    df = pd.DataFrame({
        'x': x_values,
        'y': y_values,
        'color': y_values  # Use the y-values as the color
    })

    fig = px.bar(df, x='x', y='y', color='color',
                 color_continuous_scale=['red', 'yellow', 'green'], height=300)
    fig.update_layout(autosize=False, width=2605, coloraxis_showscale=False)

    return fig


if __name__ == '__main__':
    app.run(debug=True)
