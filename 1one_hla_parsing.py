#-*-coding:utf-8
import re
import pandas as pd
from collections import OrderedDict
from json import dumps
import json
import sys

#reload(sys)
#sys.setdefaultencoding('utf-8') #한글unicode 표현 방지

################ Regular term extraction ######################
diagnosis_re = re.compile(r'Clinical Indication of study(.+)\n')

patientName1_re    = re.compile(r'Name\s+([^A-Z^\s^a-z^\d^-]+)\s+\-+', re.M)
patientName_eng_re = re.compile(r'Name(\s+)([A-Za-z]+.+)\s+-', re.M)
patientName2_re    = re.compile(r'Name\s+([^A-Z^\s^a-z^\d^-]+)\s+([^A-Z^\s^a-z^\d^-]+)\s+\-+', re.M)
patientName2_2_re  = re.compile(r'Name\s+([\S]+)\s+([A-Z]+\s+#\d+)\s+\n')
patientName3_re    = re.compile(r'Name\s+([가-힇]{1,4})\s+([가-힇]{1,4})\s+([가-힇]{1,4})\s+', re.M)
patientName4_re    = re.compile(r'Name\s+([가-힇]{1,4})\s+([가-힇]{1,4})\s+([가-힇]{1,4})\s+([가-힇]{1,4})\s+\n', re.M)

sample_no_re    = re.compile(r'Sample(\s+|\s?)No.\s+(\d+\-[A-Z]+)\s+', re.M)
sample_no2_re   = re.compile(r'Sample(\s+|\s?)No.\s+(\d+\-[A-Z]+)\s+(\d+\-[A-Z]+\(?.+\)?)\s+', re.M)
date_of_test_re = re.compile(r'Date\s+of\s+test\s+(\d+.\d+.\s?\d+)\s+', re.M) 
date_of_test2_re = re.compile(r'Date\s+of\s+test\s{18,22}(\d+.\d+.\s?\d+)\s+', re.M)
######### sex and age #########
sex_age_re = re.compile(r'((Sex\/Age)|(Age\/Sex))\s+([A-Z]+)\/(\d+)', re.M)

###############################################################################
def return_empty_clinical():
	clinical_empty_dict = {"Name":"Empty","PatientType":"Empty","diagnosis":"Empty","RelatedPatientsNo":"Empty",\
				"RelatedPatientName":"Empty","SampleNo1":"Empty", "DateOfTest1":"Empty" }

	return clinical_empty_dict

def return_empty_lab():
	lab_empty_dict = {"MRN":"Empty","PatientName":"Empty","PatientSex":"Empty","PatientDOB":"Empty",\
			"TestCode":"Empty","TestDate":"Empty","TestName":"Empty"}
	return lab_empty_dict

def return_empty_HLA():
	A_HLA_empty_dict  = {"A1_gene":"Empty","A1_allele":"Empty","A2_gene":"Empty","A2_allele":"Empty"}
	B_HLA_empty_dict  = {"B1_gene":"Empty","B1_allele":"Empty","B2_gene":"Empty","B2_allele":"Empty"}
	C_HLA_empty_dict  = {"C1_gene":"Empty","C1_allele":"Empty","C2_gene":"Empty","C2_allele":"Empty"}
	R_HLA_empty_dict = {"DR1_gene":"Empty","DR1_allele":"Empty","DR2_gene":"Empty","DR2_allele":"Empty"}
	DQ_HLA_empty_dict = {"DQ1_gene":"Empty","DQ1_allele":"Empty","DQ2_gene":"Empty","DQ2_allele":"Empty"}

	return A_HLA_empty_dict, B_HLA_empty_dict, C_HLA_empty_dict, DR_HLA_empty_dict, DQ_HLA_empty_dict

################################################################################################
def return_notMatched_HLA_A():
	A_HLA_notMatched_dict  = {"A1_gene":"Not matched","A1_allele":"Not matched","A2_gene":"Not matched","A2_allele":"Not matched"}
	return A_HLA_notMatched_dict
def return_notMatched_HLA_B():
	B_HLA_notMatched_dict  = {"B1_gene":"Not matched","B1_allele":"Not matched","B2_gene":"Not matched","B2_allele":"Not matched"}
	return B_HLA_notMatched_dict
def return_notMatched_HLA_C():
	C_HLA_notMatched_dict  = {"C1_gene":"Not matched","C1_allele":"Not matched","C2_gene":"Not matched","C2_allele":"Not matched"}
	return C_HLA_notMatched_dict
def return_notMatched_HLA_DR():
	DR_HLA_notMatched_dict = {"DR1_gene":"Not matched","DR1_allele":"Not matched","DR2_gene":"Not matched","DR2_allele":"Not matched"}
	return DR_HLA_notMatched_dict
def return_notMatched_HLA_DQ():
	DQ_HLA_notMatched_dict = {"DQ1_gene":"Not matched","DQ1_allele":"Not matched","DQ2_gene":"Not matched","DQ2_allele":"Not matched"}
	return DQ_HLA_notMatched_dict


#################################################################
def fill_mylist_fx(input_list, a, outlist):
	for items in input_list:
		outlist.append(a[items])
	return outlist
##################################################################
def get_regex_group_fx_2(input_regex):
	outlist = []
	out1 = input_regex.group(1)
	out2 = input_regex.group(2)
	outlist.append(out1)
	outlist.append(out2)
	return outlist

def get_regex_group_fx_1(input_regex):
	outlist = []
	out1 = input_regex.group(1)
	outlist.append(out1)
	return outlist

###################################################################
def HLA_A_load ( Alist ):
	result_dict = {}
	result_dict["A1_gene"]   = Alist[0]
	result_dict["A1_allele"] = Alist[1]
	result_dict["A2_gene"]   = Alist[2]
	result_dict["A2_allele"] = Alist[3]
	return result_dict

##################################################################
def HLA_B_load ( Blist ):
	result_dict = {}
	result_dict["B1_gene"]   = Blist[0]
	result_dict["B1_allele"] = Blist[1]
	result_dict["B2_gene"]   = Blist[2]
	result_dict["B2_allele"] = Blist[3]
	return result_dict

##################################################################
def HLA_C_load ( Clist ):
	result_dict = {}
	result_dict["C1_gene"]   = Clist[0]
	result_dict["C1_allele"] = Clist[1]
	result_dict["C2_gene"]   = Clist[2]
	result_dict["C2_allele"] = Clist[3]
	return result_dict

##################################################################
def HLA_DR_load ( Rlist ):
	result_dict = {}
	result_dict["DR1_gene"]   = Rlist[0]
	result_dict["DR1_allele"] = Rlist[1]
	result_dict["DR2_gene"]   = Rlist[2]
	result_dict["DR2_allele"] = Rlist[3]
	return result_dict

##################################################################
def HLA_DQ_load ( Qlist ):
	result_dict = {}
	result_dict["DQ1_gene"]   = Qlist[0]
	result_dict["DQ1_allele"] = Qlist[1]
	result_dict["DQ2_gene"]   = Qlist[2]
	result_dict["DQ2_allele"] = Qlist[3]
	return result_dict

##############222222222222222222222222222222222222222222222222222222222222222###############################
###### HLA A #######
a19_re  = re.compile(r'(HLA-[DRB1|DQB1]+)(\s+|\s?)\(?D?N?A?\)?(\s+|\s?)(\d+)\s+(\d+)', re.M)
a18_re  = re.compile(r'(HLA-[DRB1|DQB1]+)(\s+|\s?)\(?D?N?A?\)?(\s+|\s?)([DR|DQ]+(\s+|\s?)[\d+|\d?])(\s?|\s+)\*?\((.+)\)\s+([DR|DQ]+(\*|\d+))\s?\((.+)\)', re.M)
a17_re  = re.compile(r'(HLA-[DRB1|DQB1]+)(\s+|\s?)\(?D?N?A?\)?\s+(NT)', re.M)
a16_re  = re.compile(r'(HLA-[DRB1|DQB1]+)(\s+|\s?)\(?D?N?A?\)\s+(\d+)\s+(-)', re.M)
#a15_re  = re.compile(r'(HLA-[DQB1]+)(\s+|\s?)\(?D?N?A?\)\s+[DQ]+\d+\((.+)\)\s+(ND\*)\((.+)\)', re.M)

a14_re  = re.compile(r'(HLA-[C|Cw]+)(\s+|\s?)\(?D?N?A?\)?\s+([Cw]+\d+)\((.+)\)\s+(.?)\((.+)\)', re.M)
a13_re  = re.compile(r'(HLA-[C|Cw]+)(\s+|\s?)\(?D?N?A?\)?(-)\s+(-)?', re.M)


a20_re  = re.compile(r'(HL?A-[A|B|C|Cw|DRB1|DQB1]+)(\s+|\s?)\(?D?N?A?\)?\s+([A|B|C|Cw|DR|DQ]+\d+)(\s+|\s?)\((\*\d+)\)\s+([A|B|C|Cw|DR|DQ]+\d+)(\s+|\s?)\((\*\d+)\)\s+', re.M)
a21_re  = re.compile(r'(HL?A-[A|B|C|Cw|DRB1|DQB1]+)(\s+|\s?)\(?D?N?A?\)?\s+([A|B|C|Cw|DR|DQ]+\d+)(\s+|\s?)\((\*\d+:\d+g?\s?)\)\s+([A|B|C|Cw|DR|DQ]+\d+)(\s+|\s?)\((\*\d+)\)\s+', re.M)
a22_re  = re.compile(r'(HL?A-[A|B|C|Cw|DRB1|DQB1]+)(\s+|\s?)\(?D?N?A?\)?\s+([A|B|C|Cw|DR|DQ]+\d+)(\s+|\s?)\((\*\d+:\d+[a-z]?\s?)\)\s+([A|B|C|Cw|DR|DQ]+\d+)(\s+|\s?)\((\*\d+:\d+[a-z]?)\)\s+', re.M)
a23_re  = re.compile(r'(HL?A-[A|B|C|Cw|DRB1|DQB1]+)(\s+|\s?)\(?D?N?A?\)?\s+(\")\s+', re.M)
a24_re  = re.compile(r'(HL?A-[A|B|C|Cw|DRB1|DQB1]+)(\s+|\s?)\(?D?N?A?\)?\s+(([A|B|C|Cw|DR|DQ]+\d+)(\s+|\s?))\((\*\d+[a-z]?)\)\s+(-)\s+', re.M)
a25_re  = re.compile(r'(HL?A-[A|B|C|Cw|DRB1|DQB1]+)(\s+|\s?)\(?D?N?A?\)?\s+([A|B|C|Cw|DR|DQ]+\d+)(\s+|\s?)\((\*\d+)\)\s+([A|B|C|Cw|DR|DQ]+\d+)(\s+|\s?)\((\*\d+:\d+[a-z]?)\)\s+', re.M)
a26_re  = re.compile(r'(HL?A-[A|B|C|Cw|DRB1|DQB1]+)(\s+|\s?)\(?D?N?A?\)?\s+([A|B|C|Cw|DR|DQ]+\d+)(\s+|\s?)\((\*\d+:\d+g?)\)\s+(-)(\s+|\s?)\((\*\d+:\d+[a-z]?)\)\s+', re.M)
a27_re  = re.compile(r'(HL?A-[A|B|C|Cw|DRB1|DQB1]+)(\s+|\s?)\(?D?N?A?\)?\s+([A|B|C|Cw|DR|DQ]+\d+)(\s+|\s?)\((\*\d+:\d+g?)\)\s+(-)(\n+|\s+)-', re.M)
a28_re  = re.compile(r'(HL?A-[A|B|C|Cw|DRB1|DQB1]+)(\s+|\s?)\(?D?N?A?\)?\s+([A|B|C|Cw|DR|DQ]+\d+)(\s+|\s?)\((\*\d+:\d+)\)\s+([A|B|C|Cw|DR|DQ]+\d+)(\s+|\s?)\((\*\d+:\d+/\d+)\)\s+', re.M)
a29_re  = re.compile(r'(HL?A-[A|B|C|Cw|DRB1|DQB1]+)(\s+|\s?)\(?D?N?A?\)?\s+([A|B|C|Cw|DR|DQ]+\d+)(\s+|\s?)\((\*\d+:\d+)\)\s+(-)(\s+|\s?)\((\*:?\d+:\d+)\)\s+', re.M)
a30_re  = re.compile(r'(HL?A-[A|B|C|Cw|DRB1|DQB1]+)(\s+|\s?)\(?D?N?A?\)?\s+([A|B|C|Cw|DR|DQ]+\d+)\s+(-)', re.M)
a31_re  = re.compile(r'(HL?A-[A|B|C|Cw|CW|DRB1|DQB1]+)(\s+|\s?)(\(?D?N?A?\)?)(\s?|\s+)([A|B|C|Cw|CW|DR|DQ]+(\d+|\d?))\s+([A|B|C|Cw|CW|DR|DQ]+(\d+|\d?))\s+', re.M)
a32_re  = re.compile(r'(HL?A-[A|B|C|Cw|DRB1|DQB1]+)(\s+|\s?)\(?D?N?A?\)?\s+([A|B|C|Cw|DR|DQ]+)(\s+|\s?)-(\s+|\s?)\((\*\d+:\d+g?)\)\s+\((\*\d+:\d+[a-z]?)\)\s+', re.M)
a33_re  = re.compile(r'(HL?A-[A|B|C|Cw|DRB1|DQB1]+)(\s+|\s?)\(?D?N?A?\)?\s+([A|B|C|Cw|DR|DQ]+\d+)(\s+|\s?)\((\*\d+:\d+/\d+)\)\s+([A|B|C|Cw|DR|DQ]+\d+)(\s+|\s?)\((\*\d+:\d+)\)\s+', re.M)
a34_re  = re.compile(r'(HL?A-[A|B|C|Cw|DRB1|DQB1]+)(\s+|\s?)\(?D?N?A?\)?\s+([A|B|C|Cw|DR|DQ]+)(\s+|\s?)-(\s+|\s?)\((\*\d+:\d+g?)\)(\s+|\s?)-(\s+|\s?)\((\*\d+:\d+[a-z]?)\)\s+', re.M)
a35_re  = re.compile(r'(HL?A-[A|B|C|Cw|DRB1|DQB1]+)(\s+|\s?)\(?D?N?A?\)?\s+([A|B|C|Cw|c|DR|DQ]+)(\s+|\s?)(-)', re.M)
a36_re  = re.compile(r'(HL?A-[A|B|C|Cw|DRB1|DQB1]+)(\s+|\s?)\(?D?N?A?\)?\s+([A|B|C|Cw|DR|DQ]+\d+)(\s+|\s?)\((\*\d+g?\s?)\)\s+([A|B|C|Cw|DR|DQ]+\d+)(\s+|\s?)\((\*\d+g?)\)\s+', re.M)
a37_re  = re.compile(r'(HL?A-[A|B|C|Cw|DRB1|DQB1]+)(\s+|\s?)\(?D?N?A?\)?\s+([A|B|C|Cw|DR|DQ]+\d+)(\s+|\s?)\((\*\d+g?\s?)\)\s+(-)(\s+|\s?)\((\*\d+g?)\)\s+', re.M)
a38_re  = re.compile(r'(HL?A-[A|B|C|Cw|DRB1|DQB1]+)(\s+|\s?)\(?D?N?A?\)?\s+([A|B|C|Cw|DR|DQ]+\d+)(\s+|\s?)\((.+)\)\s+([A|B|C|Cw|DR|DQ]+\d+)(\s+|\s?)\((.+)\)', re.M)
a39_re  = re.compile(r'(HL?A-[A|B|C|Cw|DRB1|DQB1]+)(\s+|\s?)\(?D?N?A?\)?\s+([A|B|C|Cw|DR|DQ]+\d+)(\s+|\s?)\((.+)\)\s+([A|B|C|Cw|DR|DQ]+)(\s+|\s?)\((.+)\)', re.M)
a40_re  = re.compile(r'(HL?A-[A|B|C|Cw|DRB1|DQB1]+)(\s+|\s?)\(?D?N?A?\)?\s+([A|B|C|Cw|DR|DQ]+\d+)(\s+|\s?)\((.+)\)\s+([A|B|C|Cw|DR|DQ]+)(\s+|\s?)(\d+)', re.M)
a41_re  = re.compile(r'(HL?A-[A|B|C|Cw|DRB1|DQB1]+)(\s+|\s?)\(?D?N?A?\)?\s+([A|B|C|Cw|DR|DQ]+\d+)(\s+|\s?)\(([A|B|C|Cw|DR|DQ]?\*\d+g?\s?)\)\s+(-)(\s+|\s?)\(([A|B|C|Cw|DR|DQ]?\*\d+)\)', re.M)
a42_re  = re.compile(r'(HLA-[A|B|C|Cw|DRB1|DQB1]+)(\s+|\s?)\(?D?N?A?\)?\s+([A|B|C|Cw|DR|DQ]+\d+)(\s+|\s?)\(([A|B|C|Cw|DR|DQ]?\*\d+)\)\)?\s+(-)', re.M)
a43_re  = re.compile(r'(HL?A-[A|B|C|Cw|DRB1|DQB1]+)(\s+|\s?)\(?D?N?A?\)?\s+([A|B|C|Cw|DR|DQ]+\d+)(\s+|\s?)\((.+)\)\s+(-|Null)\s?\((.+)\)', re.M)
a44_re  = re.compile(r'(HL?A-[A|B|C|Cw|DRB1|DQB1]+)(\s+|\s?)\(?D?N?A?\)?\s+([A|B|C|Cw|DR|DQ]+\d+)(\s+|\s?)\s+([A|B|C|Cw|DR|DQ]+\d+)(\s+|\s?)\((.+)\)', re.M)
a45_re  = re.compile(r'(HL?A-[A|B|C|Cw|DRB1|DQB1]+)(\s+|\s?)\(?D?N?A?\)?\s+([A|B|C|Cw|DR|DQ]+\d+)(\s+|\s?)\((.+)\)\s+(-)', re.M)
a46_re  = re.compile(r'(HL?A-[A|B|C|Cw|DRB1|DQB1]+)(\s+|\s?)\(?D?N?A?\)?\s+([A|B|C|Cw|DR|DQ]+\*?\d+)(\s+|\s?)\((\*.+)\)\s+([A|B|C|Cw|DR|DQ]+\*?\d+)(\s+|\s?)\((\*.+)\)\s+', re.M)
a47_re  = re.compile(r'HLA-[A|B|C|Cw|DRB1|DQB1]+(\s+|\s?)\(?D?N?A?\)?\s+([A|B|C|Cw|DR|DQ]+\d+)(\s+|\s?)\((.+)\)\s+\((.+)\)\s+', re.M)
a48_re  = re.compile(r'HLA-[A|B|C|Cw|DRB1|DQB1]+(\s+|\s?)\(?D?N?A?\)?(\s+|\s?)(-)', re.M)
a50_re  = re.compile(r'(HL?A-[A|B|C|Cw|CW|DRB1|DQB1]+)(\s+|\s?)\(?D?N?A?\)?\s+(DR\d+)\s+(DR\d+)\s+',re.M)

a2NT_re = re.compile(r'(HL?A-[A|B|C|Cw|DRB1|DQB1]+)(\s+|\s?)\(?D?N?A?\)?\s+((N|n)(\s+|\s?)o?t?(\s+|\s?)(T|t)e?s?t?e?d?)\s+', re.M)

myreg_list = {a20_re:[2,4,5,7], a21_re:[2,4,5,7], a22_re:[2,4,5,7], \
		a23_re:[2,2,2,2], a24_re:[3,5,6,6], a25_re:[2,4,5,7],\
		 a26_re:[2,4,5,7], a27_re:[2,4,5,5],a2NT_re:[2,2,2,2], \
		a28_re:[2,4,5,7],a29_re:[2,4,5,7], a30_re:[2,2,3,3], a31_re:[4,4,6,6],\
		a32_re:[2,5,2,6],a33_re:[2,4,5,7],a34_re:[2,6,2,8],a35_re:[2,4,2,4],a36_re:[2,4,5,7],a37_re:[2,4,5,7],\
		a38_re:[2,4,5,7], a39_re:[2,4,5,7], a40_re:[2,4,5,5],\
		a41_re:[2,4,5,7], a42_re:[2,4,5,5], a43_re:[2,4,5,6], a44_re:[2,2,4,6],\
		a45_re:[3,4,5,5], a46_re:[2,4,5,7], a47_re:[1,3,1,4], a48_re:[1,1,1,1], a50_re:[2,2,3,3],\
		a19_re:[3,3,4,4], a18_re:[3,6,7,9], a17_re:[2,2,2,2], a16_re:[2,2,3,3],\
		a14_re:[2,3,2,5], a13_re:[2,2,2,2]}

lab_test_info_dict = {"MRN":"New","PatientName":"New","PatientSex":"New","PatientDOB":"New","TestCode":"New","TestDate":"New","TestName":"New"}
clinical_information_dict = {"Name":"New","PatientType":"New","diagnosis":"New","RelatedPatientsNo":"One","RelatedPatientName":"New","SampleNo1":"New", "DateOfTest1":"New", "TestInfo":"New", "TestRemark":"New"}

A_HLA_dict  = {"A1_gene":"New","A1_allele":"New","A2_gene":"New","A2_allele":"New"}
B_HLA_dict  = {"B1_gene":"New","B1_allele":"New","B2_gene":"New","B2_allele":"New"}
C_HLA_dict  = {"C1_gene":"New","C1_allele":"New","C2_gene":"New","C2_allele":"New"}
DR_HLA_dict = {"DR1_gene":"New","DR1_allele":"New","DR2_gene":"New","DR2_allele":"New"}
DQ_HLA_dict = {"DQ1_gene":"New","DQ1_allele":"New","DQ2_gene":"New","DQ2_allele":"New"}


test_info_re = re.compile(r'(유전자:\s?A,B,C,DQ-exon.+)')
test_info_re2 = re.compile(r'(\*\s?대상유전자\s?:.+)')

diagnosis_re = re.compile(r'Clinical Indication of study\s?:\s?(.+)')
remark_re    = re.compile(r'Remark(.+)')
remark_re2   = re.compile(r'(/*\s+참고\s+.+)')
######################################################  
#####################################################
def HLA_two_parsing(string, type):
	for myreg in myreg_list.keys():
		A = []	
		a_reg = myreg.findall(string)
		if a_reg:
			for k in range(len(a_reg)):
				if a_reg[k][0] == type:
					loc_list = myreg_list[myreg]
					A = fill_mylist_fx(loc_list, a_reg[k], A)
					return A
####################################################
############ Open Input File #########################################
input_df = pd.read_excel("/storage/home/geffa/cmi/geffa/hla/20180620_testresult_HLA_8test_supreme_83033_20180620231835_1529504315315.xlsx", \
	skiprows=[0,1],names=["ID","MRN","PatientName","Sex","DOB","TestDate","Code","TestName","TestNameDetails","DetailCode",\
				"NumericResults","BinaryResults","TestResults"])
##############################################################################################################################
state_list = ["Recipient patient with 1 donor", "Recipient patient with 2 donor","Recipient patient with 3 donor"]
####################################################################################################################################
i = 0
j = 0
input_df_size = 0 

final_list = []
unmatched  = []

Aunmatched  = []
Bunmatched  = []
Cunmatched  = []
Runmatched  = []
Qunmatched  = []

Anot  = 0
Bnot  = 0
Cnot  = 0
Rnot  = 0
Qnot  = 0
unmatched_samples = []

notmatched_test_out = open("/storage/home/geffa/cmi/geffa/hla/test_out.txt","w")

for index, row in input_df.iterrows():
	input_df_size += 1

	lab_test_info_dict = {"MRN":"New","PatientName":"New","PatientSex":"New","PatientDOB":"New","TestCode":"New","TestDate":"New","TestName":"New"}

	lab_test_info_dict["MRN"]= row['MRN']
	lab_test_info_dict["PatientName"] = row['PatientName']
	lab_test_info_dict["PatientSex"]  = row['Sex']
	lab_test_info_dict["PatientDOB"]  = row['DOB']
	lab_test_info_dict["TestDate"] = row['TestDate']
	lab_test_info_dict["TestName"] = row['TestName']
	lab_test_info_dict["TestCode"] = row['Code']
	lab_result = lab_test_info_dict
	if input_df_size < 64025:
		a = row['TestResults']
		if type(a) != float:
			data = a.replace("\n\n\n","\n")
			data = data.replace("\n\n","\n")
			result_split = data.split("\n")
			result_length = len(result_split)
			pt_name     = patientName1_re.search(data)
			pt_name_eng = patientName_eng_re.search(data)
			pt_name2    = patientName2_re.search(data)
			pt_name2_2  = patientName2_2_re.search(data)
			pt_name3    = patientName3_re.search(data)
			sampleNo    = sample_no_re.search(data)
			sampleNo2   = sample_no2_re.search(data)
			dateTest    = date_of_test_re.search(data)
			dateTest2   = date_of_test2_re.search(data)
			test_info   = test_info_re.search(data)
			test_info2  = test_info_re2.search(data)
			diagnosis   = diagnosis_re.search(data)
			remark      = remark_re.search(data)
			remark2     = remark_re2.search(data)
			if pt_name:
				j += 1
				clinical_information_dict = {"Name":"New","PatientType":"New","diagnosis":"New",\
								"RelatedPatientsNo":"One","RelatedPatientName":"New",\
								"SampleNo1":"New", "DateOfTest1":"New",\
								"TestInfo":"New", "TestRemark":"New"}

				#print (data)
				if test_info:
					ti = test_info.group(1)
					clinical_information_dict["TestInfo"] = ti
				elif test_info2:
					ti = test_info2.group(1)
					clinical_information_dict["TestInfo"] = ti
				#else:
				#	notmatched_test_out.write("%s"%data)
				#	print (data)
				if diagnosis:
					dx = diagnosis.group(1)
					clinical_information_dict["diagnosis"] = dx
				#else:
				#	notmatched_test_out.write("%s"%data)
				#	print (data)
				if remark:
					rm = remark.group(1)
					clinical_information_dict["TestRemark"] = rm
				elif remark2:
					rm = remark2.group(1)
					clinical_information_dict["TestRemark"] = rm
				#else:
				#	notmatched_test_out.write("%s"%data)
					#print (data)
				#print j
				sample_list = []

				pt_input_name = pt_name.group(1)
				if sampleNo:
					sn1 = get_regex_group_fx_2(sampleNo)[0]
				#if sampleNo2:
                                #        sn1 = get_regex_group_fx_2(sampleNo2)[0]
				#        sn2 = get_regex_group_fx_2(sampleNo2)[1]
				if dateTest:
					dt1 = get_regex_group_fx_1(dateTest)[0]
					#dt2 = get_regex_group_fx_2(dateTest)[1]
				#if dateTest2:
				#        dt1 = "NA"
				#        dt2 = get_regex_group_fx_1(dateTest2)[0]                                        

				print (pt_input_name)
				#print data
				pt_input_name  = pt_input_name.replace("-","")
				clinical_information_dict["Name"] = pt_input_name
				clinical_information_dict["PatientType"] = "Recipient patient"
				clinical_information_dict["RelatedPatientName"] = "NA"
				clinical_information_dict["SampleNo1"] = sn1
				#clinical_information_dict["SampleNo2"] = sn2
				clinical_information_dict["DateOfTest1"] = dt1
				#clinical_information_dict["DateOfTest2"] = dt2
				clinical_result = clinical_information_dict

				hla_A_result  = {}
				hla_A_result  = HLA_two_parsing(data, "HLA-A")
				hla_AnoL_result = HLA_two_parsing(data, "HA-A")
				if hla_A_result:
					A_result = HLA_A_load(hla_A_result)
				elif hla_AnoL_result:
					A_result = HLA_A_load(hla_AnoL_result)
				else:
					A_result = return_notMatched_HLA_A()
				#print A_result
			
				hla_B_result  = {}
				hla_B_result  = HLA_two_parsing(data, "HLA-B")
				if hla_B_result:
					B_result = HLA_B_load(hla_B_result)
				else:
					B_result = return_notMatched_HLA_B()
				#print B_result

				hla_C_result  = {}
				hla_C_result  = HLA_two_parsing(data, "HLA-C")
				hla_Cw_result = HLA_two_parsing(data, "HLA-Cw")
				if hla_C_result:
					C_result = HLA_C_load(hla_C_result)
				elif hla_Cw_result:
					C_result = HLA_C_load(HLA_Cw_result)
				else:
					C_result = return_notMatched_HLA_C()
				#print C_result

				hla_DR_result  = {}
				hla_DR_result  = HLA_two_parsing(data, "HLA-DRB1")
				if hla_DR_result:
					DR_result = HLA_DR_load(hla_DR_result)
				else:
					DR_result = return_notMatched_HLA_DR()
				#print DR_result

				hla_DQ_result  = {}
				hla_DQ_result  = HLA_two_parsing(data, "HLA-DQB1")
				if hla_DQ_result:
					DQ_result = HLA_DQ_load(hla_DQ_result)
				else:
					DQ_result = return_notMatched_HLA_DQ()
				#print DQ_result
				

	#		else:
	#			clinical_result = return_empty_clinical()
	#			lab_result  = return_empty_lab()
	#			A_result    = return_empty_HLA()[0]
	#			B_result    = return_empty_HLA()[1]
	#			C_result    = return_empty_HLA()[2]
	#			DR_result   = return_empty_HLA()[3]
	#			DQ_result   = return_empty_HLA()[4]
		

				sample_list.append(lab_result)
				sample_list.append(clinical_result)
				sample_list.append(A_result)
				sample_list.append(B_result)
				sample_list.append(C_result)
				sample_list.append(DR_result)
				sample_list.append(DQ_result)
				
				if A_result["A1_allele"] == "Not matched":
					Anot = Anot + 1
					unmatched_samples.append(row['MRN'])
					Aunmatched.append(row['MRN'])
					#notmatched_test_out.write("%s"%data)
					#print ("A not matched")
					#print (data)
				if B_result["B1_allele"] == "Not matched":
					Bnot = Bnot + 1
					unmatched_samples.append(row['MRN'])
					Bunmatched.append(row['MRN'])
					#notmatched_test_out.write("%s"%data)
					#print ("B not matched")
					#print (data)
				if C_result["C1_allele"] == "Not matched":
					Cnot = Cnot + 1
					unmatched_samples.append(row['MRN'])
					Cunmatched.append(row['MRN'])
					notmatched_test_out.write("%s"%data)
					print ("C not matched")
					print (data)
				if DR_result["DR1_allele"] == "Not matched":
					Rnot = Rnot + 1
					unmatched_samples.append(row['MRN'])
					Runmatched.append(row['MRN'])
					#notmatched_test_out.write("%s"%data)
					#print ("DR not matched")
					#print (data)
				if DQ_result["DQ1_allele"] == "Not matched":
					Qnot = Qnot + 1
					unmatched_samples.append(row['MRN'])
					Qunmatched.append(row['MRN'])
					#notmatched_test_out.write("%s"%data)
					#print ("DQ not matched")
					#print (data)

				final_list.append(sample_list)
	else:
        	break

print ("Read done!")

print ("Total not empty rows, n= %s"%input_df_size)
print ("Total one sample count, n= %s"%j)
print ("A not match, n= %s"%Anot)
print ("B not match, n= %s"%Bnot)
print ("C not match, n= %s"%Cnot)
print ("DR not match, n= %s"%Rnot)
print ("DQ not match, n= %s"%Qnot)


unmatched = set(unmatched_samples)
print (unmatched)


c_n_name = ['Sample','A','B','C','R','Q']
new_df = pd.DataFrame(columns = c_n_name)
new_df['Sample'] = list(unmatched )
new_df.loc[new_df['Sample'].isin(Aunmatched),'A'] =1
new_df.loc[new_df['Sample'].isin(Bunmatched),'B'] =1
new_df.loc[new_df['Sample'].isin(Cunmatched),'C'] =1
new_df.loc[new_df['Sample'].isin(Runmatched),'R'] =1
new_df.loc[new_df['Sample'].isin(Qunmatched),'Q'] =1

new_df.to_csv('OneUnmatched.csv', sep=',')


print ("Total unmatched samples, n= %s"%len(unmatched))


final_out = open("finalOutOneSample.txt", 'w')
final_out.write("Total not empty rows\tTotal one sample count\tAnotMatch\tBnotMatch\tCnotMatch\tDRnotMatch\tDQnotMatch\n")
final_out.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\n"%(input_df_size,j,Anot,Bnot,Cnot,Rnot,Qnot))


print ("Write start!")
final_index = []
numb = 0

for pt in final_list:
	numb = numb + 1
	#print numb
	#final_index.append(numb)
	result_df = pd.DataFrame()
	detail_numb = 0
	for detail in pt:
		index = [0]
		detail_numb = detail_numb + 1
		rdr = pd.DataFrame(detail, index=index)
		if detail_numb == 1:
			result_df = rdr
		else:
			result_df = pd.concat([result_df,rdr], axis = 1)
	if numb == 1:
		final_df = result_df
	else:
		final_df = pd.concat([final_df, result_df], axis=0)
#print final_df

cols = ['MRN', 'PatientDOB', 'PatientName', 'PatientSex', 'TestCode', 'TestDate', 'TestName', 'Name', 'PatientType', 'diagnosis',\
	'TestInfo','TestRemark','RelatedPatientsNo','SampleNo1','DateOfTest1',\
	 'A1_gene', 'A1_allele', 'A2_gene', 'A2_allele', \
	'B1_gene', 'B1_allele', 'B2_gene', 'B2_allele','C1_gene', 'C1_allele', 'C2_gene', 'C2_allele',\
	 'DR1_gene', 'DR1_allele', 'DR2_gene','DR2_allele', 'DQ1_gene', 'DQ1_allele', 'DQ2_gene','DQ2_allele','RelatedPatientName']

final_df = final_df[cols]


final_df.to_excel("../one_sample_only_hla_version2.xlsx", sheet_name='sheet1')

notmatched_test_out.close()
