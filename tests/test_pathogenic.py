import unittest
import os
from pathogenic import *

class TestCase(unittest.TestCase):
    def setUp(self):
        self.gnomad_path = '/media/jing/18117A5842B23232/db/gnomAD/'+\
            'release-170228'

    def tearDown(self):
        pass

    def test_get_vcfs(self):
        grange = '1:94484000-94484001'
        result = list(get_vcfs(self.gnomad_path,grange))
        self.assertEqual(result[1].split('\t')[2],'rs778234759')
        grange = '1:94461665-94461666'
        result = list(get_vcfs(self.gnomad_path,grange))
        self.assertEqual(result[0].split('\t')[2],'rs761134287')

if __name__ == '__main__':
    unittest.main()
