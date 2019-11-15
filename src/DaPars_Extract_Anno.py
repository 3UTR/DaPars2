#!/usr/bin/python

import os
import sys, getopt
import os.path


def Annotation_prepar_3UTR_extraction(gene_bed_file, gene_symbol_map_kfXref_file, output_utr_file):
    
    output_write = open(output_utr_file,'w')
    
    refseq_trapt_gene_symbol_dict = {}
    num_line = 0
    for line in open(gene_symbol_map_kfXref_file, 'r'):
        if num_line > 0:
            fields = line.strip('\n').strip('\r').split('\t')
            gene_symbol = fields[1]
            refseq_transcript_id = fields[0]
            refseq_trapt_gene_symbol_dict[refseq_transcript_id] = gene_symbol
        else:
            num_line += 1
    
    scanned_3UTR_list = []
    num_saved = 0
    for line in open(gene_bed_file,'r'):
        fields = line.strip('\n').split('\t')
        refseq_id = fields[3]
        if '_' not in fields[0]:
            
            if refseq_id not in refseq_trapt_gene_symbol_dict:
                gene_symbol = "NA"
            else:
                gene_symbol = refseq_trapt_gene_symbol_dict[refseq_id] 
            
            UTR_id = [refseq_id, gene_symbol,fields[0], fields[5]]
            UTR_id_new = '|'.join(UTR_id)
            curr_strand = fields[5]
            if curr_strand == "+":
                UTR_end = fields[2]
                gene_start = int(fields[1])
                UTR_start = str(gene_start + int(fields[-1].strip(',').split(',')[-1])+1)#1base
            elif curr_strand == "-":
                gene_start = int(fields[1])
                UTR_start = str(gene_start + 1)#1base
                UTR_end   = str(gene_start + int(fields[10].split(',')[0]))#1base, included
            
            this_UTR = fields[0]+UTR_start+UTR_end+curr_strand
            if this_UTR not in scanned_3UTR_list:
                write_line = [fields[0], UTR_start, UTR_end,UTR_id_new, '0', curr_strand]
                output_write.writelines('\t'.join(write_line) + '\n')
                scanned_3UTR_list.append(this_UTR)
                num_saved += 1
    
    
    output_write.close()   
    print "Total extracted 3' UTR: " + str(num_saved)



def Subtract_different_strand_overlap(input_gene_bed_file,output_utr_file):
    def UTRs_subtract_refine(UTRs_all):
        strand_info = UTRs_all[0].strip('\n').split('\t')[-1]
        if strand_info == '+':
            all_pos = []
            for curr_line in UTRs_all:
                left_pos = curr_line.strip('\n').split('\t')[1]
                all_pos.append(int(left_pos))
            selected_UTR_index = all_pos.index(min(all_pos))
            selected_UTR = UTRs_all[selected_UTR_index]
        else:
            all_pos = []
            for curr_line in UTRs_all:
                left_pos = curr_line.strip('\n').split('\t')[2]
                all_pos.append(int(left_pos))
            selected_UTR_index = all_pos.index(max(all_pos))
            selected_UTR = UTRs_all[selected_UTR_index]
        return selected_UTR
    temp_file = "overlap_opposite_strand_subtract.bed"
    cmd = 'subtractBed -a %s -b %s -S > %s' % (input_gene_bed_file, input_gene_bed_file, temp_file)
    os.system(cmd)
    
    read_subtract_result_dict = {}
    for line in open(temp_file,'r'):
        transcript_id = line.split('\t')[3].split('|')[0]
        if transcript_id not in read_subtract_result_dict:
            read_subtract_result_dict[transcript_id] = []
        read_subtract_result_dict[transcript_id].append(line)
    
    output_utr_write = open(output_utr_file,'w')
    for curr_trans_id in read_subtract_result_dict:
        curr_3UTRs = read_subtract_result_dict[curr_trans_id]
        num_3UTRs = len(curr_3UTRs)
        if num_3UTRs == 1:
            output_utr_write.writelines(curr_3UTRs[0])
        else:
            selected_UTR = UTRs_subtract_refine(curr_3UTRs)
            output_utr_write.writelines(selected_UTR)
    output_utr_write.close()
    
    try:
        os.remove(temp_file)
    except OSError:
        pass

def Extract_Anno_main(argv):
    gene_bed_file = ''
    gene_symbol_annotation_file = ''
    output_extract_file = 'temp_anno_extracted.bed'
    output_final_extract_file = ''
    
    try:
        opts, args = getopt.getopt(argv,"hb:s:o:",["bed=","symbol=","ofile"])
    except getopt.GetoptError:
        print 'python DaPars_Extract_Anno.py -b <gene_bed_file> -s <gene_symbol_map> -o <output_file>'
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print 'python DaPars_Extract_Anno.py -b <gene_bed_file> -s <gene_symbol_map> -o <output_file>'
            sys.exit()
        elif opt in ("-b", "--bed"):
            gene_bed_file = arg
        elif opt in ("-s", "--symbol"):
            gene_symbol_annotation_file = arg
        elif opt in ("-o", "--ofile"):
            output_final_extract_file = arg
    
    if gene_bed_file=='':
        print >> sys.stderr, "Error: No gene bed file!"
        exit(1)
    if gene_symbol_annotation_file=='':
        print >> sys.stderr, "Error: No gene symbol file!"
        exit(1)
    
    if output_final_extract_file=='':
        print >> sys.stderr, "Error: No output file!"
        exit(1)
    
    print "Generating regions ..."
    Annotation_prepar_3UTR_extraction(gene_bed_file, gene_symbol_annotation_file,output_extract_file)
    Subtract_different_strand_overlap(output_extract_file,output_final_extract_file)
    
    try:
        os.remove(output_extract_file)
    except OSError:
        pass

    
    print "Finished"

if __name__ == '__main__':
    Extract_Anno_main(sys.argv[1:])
    
