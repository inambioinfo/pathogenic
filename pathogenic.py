'''
search for all ABCA4 variants from gnomAD, to find pathogenic variants,
especially the non-coding ones. then do some population analysis using
all variants as a background distribution

python3
'''
import sys
import subprocess
sys.path.append('../gnomad')
import os

def get_vcfs(gnomad_path, grange):
    '''
    return a generator, first exome vcf, second gnome vcf
    '''
    # first get exomes
    exome_file = os.path.join(
            gnomad_path,
            'exomes',
            'vcf',
            'gnomad.exomes.r2.0.1.sites.vcf.gz'
            )
    result_exome = subprocess.check_output([
        'tabix',
        exome_file,
        grange
        ]).decode('utf8')
    
    # now genome
    chrom = grange.split(':')[0]
    genome_file = os.path.join(
            gnomad_path,
            'genomes',
            'vcf',
            'gnomad.genomes.r2.0.1.sites.{}.vcf.gz'.format(chrom),
            )
    return (
            subprocess.check_output([
                'tabix',
                f,
                grange
                ]).decode('utf8')
            for f in (exome_file,genome_file)
            )


if __name__ == '__main__':
    gnomad_path = '/media/jing/18117A5842B23232/db/gnomAD/release-170228'
    grange = '1:94458393-94586688'
    clinvar_file = '/media/jing/18117A5842B23232/db/clinvar/clinvar_20171029.vcf.gz'
