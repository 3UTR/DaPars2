import pandas as pd
import portion as P
import argparse

def main(args):
    colnames = ['chr','start', 'end', 'id', 'nb', 'strand']
    df = pd.read_csv(args.input, sep = '\t', names = colnames, header=None)
    df[['transcript','gene','chr','strand']] = df['id'].str.rsplit('|', expand=True)
    df['interval']=df.apply(lambda x:[x['start'],x['end']],axis=1)
    d = dict(zip(df['id'],df['interval']))

    remove = []

    keys = list(d.keys())
    for i in range(len(d)-1):
        i1 = d[keys[i]]
        i2 = d[keys[i+1]]
        overlap = P.closed(i1[0],i1[1]).overlaps(P.closed(i2[0],i2[1]))
        if overlap:
            transcript1,gene1,chrom1,strand1 = keys[i].split('|')
            transcript2,gene2,chrom2,strand2 = keys[i+1].split('|')
            if gene1 == gene2:
                continue
            else:
                remove.append(gene1)
                remove.append(gene2)

    df = df[~df['gene'].isin(remove)]
    df = df[colnames]
    df.to_csv(args.output, sep = '\t', index = False, header = False)


def parse_args():
    parser = argparse.ArgumentParser(
        prog='Dapars_Filter_Anno.py', 
        usage='python Dapars_Filter_Anno.py [options]',
        description='Gencodev36_3UTR_annotation'
    )
    parser.add_argument(
        '--input', type=str, default='Gencodev36_3UTR_annotation.bed',
        help=''
    )
    parser.add_argument(
        '--output', type=str, default='Gencodev36_3UTR_annotation_overlapping_filtered.bed',
        help=''
    )
    args = parser.parse_args()
    return args

if __name__ == '__main__':
    args = parse_args()
    main(args)
