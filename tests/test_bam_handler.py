import unittest
from unittest.mock import patch, MagicMock

from modiview.bam_handler import BAMFileHandler


class TestBAMFileHandler(unittest.TestCase):
    def setUp(self):
        self.filepath = 'filepath'

    @patch('pysam.AlignmentFile')
    def test_init(self, mock_alignment_file):
        BAMFileHandler(self.filepath)
        mock_alignment_file.assert_called_once_with(self.filepath, "rb", check_sq=False)

    @patch('pysam.AlignmentFile')
    def test_get_read(self, mock_alignment_file):
        mock_read = MagicMock()
        mock_alignment_file.return_value.fetch.return_value = [mock_read]
        bam_handler = BAMFileHandler(self.filepath)
        self.assertEqual(bam_handler.get_read(1), mock_read)

    @patch('pysam.AlignmentFile')
    def test_get_read_id(self, mock_alignment_file):
        mock_read = MagicMock()
        mock_read.query_name = 'read1'
        mock_alignment_file.return_value.fetch.return_value = [mock_read]
        bam_handler = BAMFileHandler(self.filepath)
        self.assertEqual(bam_handler.get_read_id(1), 'read1')

    @patch('pysam.AlignmentFile')
    def test_get_read_count(self, mock_alignment_file):
        mock_read = MagicMock()
        mock_alignment_file.return_value.fetch.return_value = [mock_read, mock_read]
        bam_handler = BAMFileHandler(self.filepath)
        self.assertEqual(bam_handler.get_read_count(), 2)

    @patch('pysam.AlignmentFile')
    def test_get_fasta_sequence(self, mock_alignment_file):
        mock_read = MagicMock()
        mock_read.query_name = 'read1'
        mock_read.query_sequence = 'ATGC'
        mock_alignment_file.return_value.fetch.return_value = [mock_read]
        bam_handler = BAMFileHandler(self.filepath)
        expected_fasta = ">read1\nATGC"
        self.assertEqual(bam_handler.get_fasta_sequence(1), expected_fasta)