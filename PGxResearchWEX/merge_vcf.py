
import os
import glob

files = glob.glob("BRONJ_*.vcf")

for f in files:
    cmd = "bgzip {}; tabix -p vcf {}.gz".format(f, f)
    print(cmd)
    os.system(cmd)

cmd = "vcf-merge BRONJ_*.vcf.gz > raw_merge_BRONZ.vcf"
print(cmd)
os.system(cmd)
