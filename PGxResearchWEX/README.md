# PGxResearchWEX
Pharmacogenomics Research in Case vs Control Setting using WEX

  << 순서 >>
1) VCF merge 하기 

1-1) 만약 index 가 필요하면 tabix 로 indexing 한다. 
  
    index_tabix.py
    
    
1-2) 이때 header 가 없는 경우 만들어서 넣어주고, GATK 를 이용한 variant filter 를 적용한다. 
 
    vcf_header.py (header file: header.vcf) 

 

1-3) 개별 vcf.gz 이라면 case, control 각각의 group 으로 묶는다.
  
    raw_merge_vcf.py 
    

1-4) 이중 snp 과 indel 을 각각 다른 파일로 나눈다. 

    extract_raw_variant.sh
    
   
1-5) case group 과 control group 을 하나의 파일로 묶는다. 

    #bgzip 으로 압축 후 tabix로 index 먼저 하고 ==>  bgzip Tacro_merged.vcf ; tabix -p vcf Tacro_merged.vcf.gz
    merge-vcf.sh

     
2) snpEFF annotation 을 수행한다. 

2-1) snpEFF annotation

    annotate.sh

2-2) snpSIFT 를 이용하여 case, control 군을 표시한 테이블로 만든다. 

    case_control.py
    
2-3) SIFT score annotation from the SIFT database합친 파일에 대해 SIFT data base 를 이용하여 SIFT score 를 annotation 한다. 

    SIFT_vcf.py
 
2-4) snpsift 의 case control study 결과를 csv 파일로 parsing 한다.

    parse_case_control.py
    
    
    
3) gene score계산하기

3-1) gene score 를 계산한다.

    gene_score.py
    
3-2) wilcoxon test 로 case vs control 차이나는 gene 을 확인한다. 

    gene_score_wilcox.py
    
3-3) 결과 vcf 를 text file 로 변환저장한다. 

    이때 먼저 header 를 만든다. 
    
        df = open("/home/genome/Tacro/merge_Tacro_PEP.vcf","r")
        header = open("/home/genome/Tacro/Tacro_PEP_header.txt","w")
        for line in df.readlines():
            if line[0] =="#" and line[1] =="C":
                header.write(line)
        header.close()
        
    parse_vcf_to_txt.py 

    
    
    
    
