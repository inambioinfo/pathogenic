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
        grange = '1:94476950-94476951'
        result = list(get_vcfs(self.gnomad_path,self.fasta_ref,grange))
        vcf_string = '##fileformat=VCFv'
        self.assertTrue(result[1].split('\n')[0].startswith(vcf_string))
    def test_parse_vcfs(self):
        grange = '1:94476900-94477000'
        data = get_vcfs(self.gnomad_path,self.fasta_ref,grange)
        data = parse_vcfs(data)
        self.assertEqual(data.loc['1-94476951-A-G']['ID'],'rs1800728')
        self.assertEqual(
                data.loc['1-94476951-A-G']['CLIN_SIG'],
                'uncertain_significance&likely_pathogenic'
                )

if __name__ == '__main__':
    unittest.main()
