import pysam


class BAMFileHandler:
    def __init__(self, filepath):
        self.file = pysam.AlignmentFile(filepath, "rb", check_sq=False)
        self.reads = list(self.file.fetch(until_eof=True))

    def get_read(self, read_number):
        return self.reads[read_number - 1]

    def get_read_id(self, read_number):
        return self.get_read(read_number).query_name

    def get_read_count(self):
        return len(self.reads)

    # Workaround because AlignmentChart takes only FASTA/ClustalW sequences
    def get_fasta_sequence(self, read_number):
        read_id = self.get_read_id(read_number)
        read_sequence = self.get_read(read_number).query_sequence

        return f">{read_id}\n{read_sequence}"
