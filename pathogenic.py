'''
search for all ABCA4 variants from gnomAD, to find pathogenic variants,
especially the non-coding ones. then do some population analysis using
all variants as a background distribution

python3
'''
import sys
import re
import subprocess
sys.path.append('../gnomad')
import os
from collections import defaultdict
import pandas as pd

def split_iter(string):
    '''
    iter through string separated by \n
    '''
    return (x.group(0) for x in re.finditer(r"[^\n]+",string))

def get_vcfs(gnomad_path, fasta_ref, grange):
    '''
    return a generator, first exome vcf, second gnome vcf
    '''
    # define the pipelines
    def pl(f):
        ps1 = subprocess.Popen((
            'tabix',
            '-h',
            f,
            grange
            ),
            stdout = subprocess.PIPE
            )
        ps2 = subprocess.Popen((
            'bcftools',
            'norm',
            '-Ou',
            '-m',
            '-any'),
            stdin = ps1.stdout,
            stdout = subprocess.PIPE
            )
        xome = subprocess.check_output((
            'bcftools',
            'norm',
            '-Ov',
            '-f',
            fasta_ref,
            ),
            stdin = ps2.stdout
            ).decode('utf8')
        ps1.wait()
        ps2.wait()
        return xome
    # get files
    chrom = grange.split(':')[0]
    files = {
            'exomes':os.path.join(
                gnomad_path,
                'exomes',
                'vcf',
                'gnomad.exomes.r2.0.1.sites.vcf.gz'
                ),
            'genomes':os.path.join(
                gnomad_path,
                'genomes',
                'vcf',
                'gnomad.genomes.r2.0.1.sites.{}.vcf.gz'.format(chrom),
                )
    }
    
    return (pl(f) for f in files.values())

def parse_vcfs(vcf_generator):
    '''
    return a pd.DataFrame, with variant ids on the index column
    columns:
    id(dbsnp id), filter, consequence(on main transcript),  ac, af, an, hom
    afr/amr/eas/fin/nfe/oth/sas (ac,an,af,hom), 
    male_ac/an/hom, female_ac/an/hom, 

    ignore asj for the time being as it is missing in exome
    '''
    vcf_list = defaultdict(dict)
    for gg in vcf_generator:
        header = []
        for row in split_iter(gg):
            # skip comments and get header
            if row.startswith('##'): continue
            row = row.split('\t')
            if row[0].startswith('#'):
                header = row
                print(header)
                continue
            v_id = '-'.join([row[0],row[1],row[3],row[4]])
            for field in ('filter','id'):
                vcf_list[v_id][field] = row[header.index(field.upper())]
    print(vcf_list)


if __name__ == '__main__':
    gnomad_path = '/media/jing/18117A5842B23232/db/gnomAD/release-170228'
    grange = '1:94458393-94586688'
    clinvar_file = '/media/jing/18117A5842B23232/db/clinvar/clinvar_20171029.vcf.gz'
    fasta_ref = '/media/jing/18117A5842B23232/db/human_g1k_v37.fasta'
