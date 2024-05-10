from dash import Dash, html, dcc, callback, Output, Input
import pysam
import dash_bio as dashbio

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
    html.Div(id='slider-output-container'),
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
    ),
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
    read_id = get_read_id(value)
    return get_fasta_sequence(read_id)


if __name__ == '__main__':
    app.run(debug=True)
