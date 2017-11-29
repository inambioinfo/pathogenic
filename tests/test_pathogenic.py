import unittest
import os
from pathogenic import *
import subprocess

class TestCase(unittest.TestCase):
    def setUp(self):
        self.gnomad_path = 'tests/data/gnomad'
        fasta_gz = 'tests/data/chrom1.fasta.gz'
        self.fasta_ref = 'tests/data/chrom1.fasta'
        if not os.path.isfile(self.fasta_ref):
            with open('tests/data/chrom1.fasta','w') as outf:
                subprocess.run((
                        'gunzip',
                        '-c',
                        fasta_gz,
                        ),stdout = outf)
                


    def tearDown(self):
        pass

    def test_get_vcfs(self):
        grange = '1:94484000-94484001'
        result = list(get_vcfs(self.gnomad_path,self.fasta_ref,grange))
        self.assertEqual(result[1].split('\t')[2],'rs778234759')
        grange = '1:94461665-94461666'
        result = list(get_vcfs(self.gnomad_path,self.fasta_ref,grange))
        self.assertEqual(result[0].split('\t')[2],'rs761134287')
    def test_parse_vcfs(self):
        grange = '1:94461636-94461666'
        vcfs = get_vcfs(self.gnomad_path,self.fasta_ref,grange)
        result = parse_vcfs(vcfs)
        print(result)

if __name__ == '__main__':
    unittest.main()
