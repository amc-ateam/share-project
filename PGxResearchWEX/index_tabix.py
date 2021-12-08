import os
import glob

files = glob.glob("BRONJ_*.vcf")

for f in files:
    cmd = "bgzip {}; tabix -p vcf {}.gz".format(f, f)
    print(cmd)
    os.system(cmd)
