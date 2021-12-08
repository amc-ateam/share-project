
import gzip
import re
import sys
import vcf
import MySQLdb
import pandas as pd
import scipy.stats
from scipy import stats
import numpy as np
from tqdm import tqdm

# Using Factory Pattern
# https://www.packtpub.com/books/content/python-design-patterns-depth-factory-pattern

class VCFConnector:

    def __init__(self, filepath):
        self.data = []

        # PyVCF: https://pyvcf.readthedocs.io/en/latest/INTRO.html
        vcf_reader = vcf.Reader( open(filepath, 'r') )

        for record in vcf_reader:
            if record.is_snp and ( (not record.FILTER) or record.FILTER=='' or record.FILTER is None or record.FILTER=='PASS' ):
                self.data.append( {'chrom':record.CHROM, 'pos':record.POS, 'ref':record.REF, 'alt':record.ALT[0]} )

    @property
    def parsed_data(self):
        return self.data

def connection_factory(filepath):
    if filepath.endswith('.vcf') or filepath.endswith('vcf.gz'):
        connector = VCFConnector
    else:
        raise ValueError('Cannot connect to {}'.forpath(filepath))

    return connector(filepath)

def connect_to(filepath):
    factory = None
    try:
       factory = connection_factory(filepath)
    except ValueError as ve:
       print ve

    return factory

import glob
import re

files = glob.glob("header/pass_*.vcf.gz")

con = MySQLdb.connect('localhost', 'root', 'snubi$manse#', 'xentinel')
cur = con.cursor(MySQLdb.cursors.DictCursor)

sift_scores_all = pd.DataFrame([], columns = ['chrom', 'pos', 'ref', 'alt', 'symbol', 'sift_score', 'sample_id'])
gene_scores_all = pd.DataFrame([], columns=['symbol', 'gene_score', 'sample_id'])

for filename in tqdm(files):
    # header/pass_filtered_BRONJ_00015.vcf.gz
    m = re.search(r"header/pass_filtered_(\S+)\.vcf\.gz", filename)
    sample_id = m.group(1)

    vcf_factory = connect_to(filename)
    person_data = vcf_factory.parsed_data

    #print "Sift Score Mapping.."

    sift_scores = []
    for row in person_data:
        chrom = row['chrom']
        pos = row['pos']
        ref = row['ref']
        alt = row['alt']

        sql = "SELECT symbol, ref, alt, score FROM sift \
               WHERE chrom='%s' AND pos='%s' AND type != 'Reference' AND ref='%s' AND alt='%s'" % ( chrom, pos, ref, alt )
        cur.execute(sql)
        rows = cur.fetchall()
        for r in rows:
            sift_scores.append( [chrom, pos, ref, alt, r['symbol'], r['score'] ] )

    sift_scores = pd.DataFrame(sift_scores)
    sift_scores.columns = ['chrom', 'pos', 'ref', 'alt', 'symbol', 'sift_score']
    sift_scores['sift_score'].replace(0, 0.00000001, inplace=True)
    #sift_scores = sift_scores[ sift_scores['sift_score'] < 0.7 ]
    sift_scores['sample_id'] = sample_id

    sift_scores_all = pd.concat([sift_scores_all, sift_scores], ignore_index=True)
    #print "Insert Sift Score.."
    #cur.executemany("INSERT INTO sift_scores(chrom, pos, ref, alt, symbol, sift_score, auth_user_id) VALUES (%s, %s, %s, %s, %s, %s, %s)", sift_scores.values.tolist() )
    #con.commit()

    # Calc Variant Scores
    variant_scores = sift_scores.groupby( ['chrom', 'pos', 'ref', 'alt', 'symbol'], as_index=False  )['sift_score'].min()

    # Calc Gene Scores
    #print "Calc Gene Score.."
    gene_scores = variant_scores.groupby( ['symbol'] )[ 'sift_score' ].apply( scipy.stats.mstats.gmean )
    gene_scores = gene_scores.reset_index()
    gene_scores.columns = ['symbol', 'gene_score']
    gene_scores['sample_id'] = sample_id

    gene_scores_all = pd.concat([gene_scores_all, gene_scores], ignore_index=True)

sift_scores_all.to_pickle("BRONJ_sift_scores.pkl")
gene_scores_all.to_pickle("BRONJ_gene_scores.pkl")
