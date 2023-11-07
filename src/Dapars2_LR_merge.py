import glob
import pandas as pd
import argparse
from statsmodels.sandbox.stats.multicomp import multipletests


pd.set_option('display.float_format', '{:.2e}'.format)

def reject(df,x):
	if float(df['p-value']) <0.05:
		if abs(float(df['DeltaPDUI'])) > x:
			return 'No'
		else:
			return 'Yes'
	elif float(df['p-value']) >=0.05:
		return 'Yes'
	else:
		return 'Yes'

def main(args):
	infile = glob.glob("./*_chr*/*.txt")

	first = True

	for f in infile:
		df = pd.read_csv(f, sep = '\t')
		if not first:
			df_final = pd.concat([df_final, df], axis=0)
		if first:
			df_final = df
			first = False

	df_final = df_final[~df_final['p-value'].isna()]

	df_final['DeltaPDUI'] = df_final[df_final.columns[-3]] - df_final[df_final.columns[-2]]

	df_final['Gene'] = df_final['Gene'].str.split('|', expand=True)[1]

	df_final = df_final.sort_values('DeltaPDUI',key=abs).drop_duplicates('Gene', keep='last')

	df_final['p-value'] = multipletests(df_final['p-value'], method='fdr_bh')[1]

	df_final['Reject'] =  df_final.apply(reject,x=args.delta,axis = 1)

	df_final.to_csv(args.out, sep = '\t', index = False)

def parse_args():
	parser = argparse.ArgumentParser(
		prog='merge_Dapars.py', 
		usage='python3 merge_Dapars.py -o OutFile --delta PDUI_significant'
	)
	parser.add_argument(
		'--out','-o', type=str, default='Dapars2_all_chr.txt'),
	parser.add_argument(
		'--delta', type=float, default=0.1)

	
	args = parser.parse_args()
	return args

if __name__ == '__main__':
	args = parse_args()
	main(args)


