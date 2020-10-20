import pandas as pd
import numpy as np

head = ['study','case','biospecimen','sample','aliquot','analyte','cell_image','read_group','mass_cytometry_assay','mass_cytometry_image','read_group_qc','submitted_unaligned_read','submitted_somatic_mutation','submitted_copy_number']
nodes = ["studies","cases","biospecimens","samples","aliquots","analytes","mass_cytometry_assays","cell_images","read_groups","mass_cytometry_images","read_group_qcs","submitted_unaligned_reads","submitted_somatic_mutations","submitted_copy_numbers"]
s_nodes = ["immunoassays","quantification_assays"]
head_s_nodes = ["immunoassay","quantification_assay"]

def hasNextNode(data):
	child = list((set(data.keys()) & set(nodes) )  | (set(data.keys()) & set(head) ) )
	return ( len(child) > 0 )

def hasSNode(data):
	child = list((set(data.keys()) & set(s_nodes) ) | (set(data.keys()) & set(head_s_nodes) ) )
	return ( len(child) > 0)



### start json to df conversion ########
def jsonToDF(data):
	child = list(  (set(data.keys()) & set(nodes) )  | (set(data.keys()) & set(head) ) )
	#hchild = list( set(data.keys()) & set(head) )
	cols = list( set(child) ^ set(data.keys())  )
	cols = [n for n in cols if n not in s_nodes]
	cols = [n for n in cols if n not in head_s_nodes]
	df = pd.DataFrame()
	if hasNextNode(data):
		for childname in child:
			for i in range(len(data[childname])):
				#combine df together!!!! #df = jsonToDF(data[childname][i])'
				tmp = jsonToDF(data[childname][i])
				df = df.append(tmp, ignore_index=True)
	di = len(df.index)
	if di == 0:
		df = pd.DataFrame(index=range(1))
	for col in cols:
		df[col] = data[col]
	if hasSNode(data):
		#### handle getting analyte info here as columns within df
		child = list(  (set(data.keys()) & set(s_nodes) )  | (set(data.keys()) & set(head_s_nodes) ) )
		for childname in child:
			if childname in ['quantification_assay','quantification_assays'] :
				for qa in range(len(data[childname])):
					for qcol in data[childname][qa].keys():
						df[qcol] = data[childname][qa][qcol]
			if childname in ['immunoassay','immunoassays'] :
				for qa in range(len(data[childname])):
					for qcol in data[childname][qa].keys():
						df[qcol] = data[childname][qa][qcol]
	return df
