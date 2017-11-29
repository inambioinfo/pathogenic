'''
search for all ABCA4 variants from gnomAD, to find pathogenic variants,
especially the non-coding ones. then do some population analysis using
all variants as a background distribution

python3
'''
import sys
import re
import subprocess
import itertools
sys.path.append('../gnomad')
import os
from collections import defaultdict,namedtuple
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

def parse_vcfs(vcf_generator, PASS=True):
    '''
    return a pd.DataFrame, with variant ids on the index column
    columns:
    id(dbsnp id), filter, consequence(on main transcript),  ac, af, an, hom
    afr/amr/eas/fin/nfe/oth/sas (ac,an,af,hom), 
    male_ac/an/hom, female_ac/an/hom, 

    ignore asj for the time being as it is missing in exome

    if PASS is True, only look at variants that pass the filter.
    '''
    vcf_list = defaultdict(dict)

    # get info fields to extract
    info_fields = (
            'AC',
            'AN',
            'AF',
            'Hom',
            )
    pop = (
            'AFR',
            'AMR',
            'EAS',
            'FIN',
            'NFE',
            'OTH',
            'SAS',
            'Male',
            'Female',
            )
    pop_info = ('_'.join(i) for i in itertools.product(info_fields,pop))
    info_fields = info_fields + tuple(pop_info)
    for gg in vcf_generator:
        for row in split_iter(gg):
            # skip comments and get row header and CSQ header
            if row.startswith('##'):
                if row.startswith('##INFO=<ID=CSQ'):
                    header = row.split('"')[1].split('Format: ')[1].split('|')
                    CSQ = namedtuple('CSQ',header)
                continue
            row = row.split('\t')
            if row[0].startswith('#'):
                header = row
                # convert header into a namedtuple class
                header[0] = header[0][1:]
                Record = namedtuple('Record',header)
                continue
            record = Record(*row)
            # not interested in variants that don't pass the filter
            if PASS and record.FILTER != 'PASS': continue

            # get variant_id
            v_id = '-'.join((
                record.CHROM,
                record.POS,
                record.REF,
                record.ALT,
                ))

            # get FILTER and ID. replace ID's '.' with None
            vcf_list[v_id]['FILTER'] = record.FILTER
            vcf_list[v_id]['ID'] = record.ID != '.' or None

            # deal with pop info
            info = record.INFO.split(';')
            info_dict = {}
            for i in info:
                if '=' not in i: continue
                ii = i.split('=')
                info_dict[ii[0]] = ii[1]
            for pop in info_fields:
                if pop in info_dict:
                    if 'AF' in pop:
                        vcf_list[v_id][pop] = float(info_dict[pop])
                    else:
                        vcf_list[v_id][pop] = int(info_dict[pop])
                else:
                    vcf_list[v_id][pop] = None

            # deal with CSQ. only parse the first annotation.
            first_csq = info_dict['CSQ'].split(',')[0].split('|')
            csq = CSQ(*first_csq)
            csq_fields = (
                    'Consequence',
                    'IMPACT',
                    'SYMBOL',
                    'Gene',
                    'Feature',
                    'HGVSc',
                    'HGVSp',
                    'EXON',
                    'INTRON',
                    'CLIN_SIG',
                    )
            for field in csq_fields:
                vcf_list[v_id][field] = getattr(csq, field)
    print(pd.DataFrame(vcf_list).T)


if __name__ == '__main__':
    gnomad_path = '/media/jing/18117A5842B23232/db/gnomAD/release-170228'
    grange = '1:94458393-94586688'
    clinvar_file = '/media/jing/18117A5842B23232/db/clinvar/clinvar_20171029.vcf.gz'
    fasta_ref = '/media/jing/18117A5842B23232/db/human_g1k_v37.fasta'
    
    vcfs = get_vcfs(gnomad_path,fasta_ref,grange)
    result = parse_vcfs(vcfs)
