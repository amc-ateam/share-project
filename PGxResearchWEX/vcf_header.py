
import glob
import gzip
import os
from tqdm import tqdm

with open("header.vcf", "r") as f:
    header = f.read()

for filename in tqdm(glob.glob("*.vcf.gz")):
    sample_id = filename.split('.')[0]
    with gzip.open(filename, 'rt') as f:
        content = f.read()

    vcf_header = '#CHROM        POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  {}\n'.format(sample_id)

    content = header + vcf_header + content

    with open("header/{}.vcf".format(sample_id), "w") as fout:
        fout.write(content)

    cmd = '''
        java -Xmx8G -jar $GATK \
        -R $REFERENCE_GENOME \
        -T SelectVariants \
        --variant header/{}.vcf \
        -o header/snp_{}.vcf \
        -selectType SNP
    '''.format(sample_id, sample_id)
    res = os.system(cmd)

    cmd = '''
        java -Xmx8G -jar $GATK \
        -R $REFERENCE_GENOME \
        -T SelectVariants \
        --variant header/{}.vcf \
        -o header/indel_{}.vcf \
        -selectType INDEL
    '''.format(sample_id, sample_id)
    res = os.system(cmd)

    cmd = '''
          java -Xmx8G -jar $GATK \
          -R $REFERENCE_GENOME \
          -T VariantFiltration \
          --variant header/snp_{}.vcf \
          -o header/filtered_snp_{}.vcf \
          --filterExpression "DP < 10" --filterName "DP10" \
          --filterExpression "MQ < 40.0" --filterName "MQ40" \
          --filterExpression "FS > 60.0" --filterName "FS60"
    '''.format(sample_id, sample_id)
    res = os.system(cmd)

    cmd = '''
          java -Xmx8G -jar $GATK \
          -R $REFERENCE_GENOME \
          -T VariantFiltration \
          --variant header/indel_{}.vcf \
          -o header/filtered_indel_{}.vcf \
          --filterExpression "DP < 10" --filterName "DP10" \
          --filterExpression "MQ < 40.0" --filterName "MQ40" \
          --filterExpression "FS > 200.0" --filterName "FS200"
    '''.format(sample_id, sample_id)
    res = os.system(cmd)

    cmd = '''
        java -Xmx8G -jar $GATK \
       -T CombineVariants \
       -R $REFERENCE_GENOME \
       --variant:snp header/filtered_snp_{}.vcf \
       --variant:indel header/filtered_indel_{}.vcf \
       -o header/filtered_{}.vcf \
       --genotypemergeoption PRIORITIZE \
       -priority snp,indel
    '''.format(sample_id, sample_id, sample_id)
    res = os.system(cmd)

    cmd = '''
          java -Xmx8G -jar $GATK \
          -R $REFERENCE_GENOME \
          -T SelectVariants \
          --variant header/filtered_{}.vcf \
          -select 'vc.isNotFiltered()' \
          -o header/pass_filtered_{}.vcf
    '''.format(sample_id, sample_id)
    res = os.system(cmd)

    os.system("bgzip header/pass_filtered_{}.vcf".format(sample_id))
