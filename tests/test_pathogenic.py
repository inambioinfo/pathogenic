import unittest
import os
from pathogenic import *

class TestCase(unittest.TestCase):
    def setUp(self):
        self.gnomad_path = '/media/jing/18117A5842B23232/db/gnomAD/'+\
            'release-170228'
        self.fasta_ref = '/media/jing/18117A5842B23232/db/human_g1k_v37.fasta'

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
