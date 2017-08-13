import os
import sys
import datetime
import threading
import numpy as np
from multiprocessing import Pool
from bisect import bisect
import time
import multiprocessing


def time_now():
    curr_time = datetime.datetime.now()
    return curr_time.strftime("%c")


def Convert_wig_into_bp_coverage(extracted_coverage, extracted_3UTR_region, strand_info):
    bp_coverage = np.zeros(extracted_3UTR_region[-1] - extracted_3UTR_region[0])
    relative_start = extracted_3UTR_region[0]
    for i in range(len(extracted_coverage)):
        curr_region_start = extracted_3UTR_region[i] - relative_start
        curr_region_end = extracted_3UTR_region[i+1] - relative_start
        bp_coverage[curr_region_start:curr_region_end] = extracted_coverage[i]
    if strand_info == '-':
        bp_coverage = bp_coverage[::-1]
    return bp_coverage


def parse_cfgfile(cfg_file):
    Group1_Tophat_aligned_file = ''
    Group2_Tophat_aligned_file = ''
    output_directory = ''
    Annotated_3UTR_file = ''
    Output_result_file = ''
    Num_least_in_group1_local = ''
    Num_least_in_group2_local = ''
    Coverage_cutoff_local = ''
    Num_threads = 1
    sequencing_depth_file = ''

    for line in open(cfg_file, 'r'):
        if line[0] == '\n' or line[0] == '#':
            comments = line
        else:
            line = line.rstrip()
            command = line.split('=')
            if command[0] == 'Group1_Tophat_aligned_Wig':
                Group1_Tophat_aligned_file = command[1].split(',')
            if command[0] == 'Group2_Tophat_aligned_Wig':
                Group2_Tophat_aligned_file = command[1].split(',')
            if command[0] == 'Output_directory':
                output_directory = command[1]
                if output_directory[-1] != '/':
                    output_directory += '/'
            if command[0] == 'Annotated_3UTR':
                Annotated_3UTR_file = command[1]
            if command[0] == 'Output_result_file':
                Output_result_file = command[1]

            # Parameters
            if command[0] == 'Num_least_in_group1':
                Num_least_in_group1_local = command[1]
            if command[0] == 'Num_least_in_group2':
                Num_least_in_group2_local = command[1]
            if command[0] == 'Coverage_cutoff':
                Coverage_cutoff_local = command[1]
            if command[0] == 'Num_Threads':
                Num_threads = int(command[1])
            if command[0] == 'sequencing_depth_file':
                sequencing_depth_file = command[1]

    if Group1_Tophat_aligned_file == '':
        print >> sys.stderr, "No Tophat aligned BAM file for group 1!"
        exit(1)
    if Group2_Tophat_aligned_file == '':
        print >> sys.stderr, "No Tophat aligned BAM file for group 2!"
        exit(1)
    if output_directory == '':
        print >> sys.stderr, "No output directory!"
        exit(1)
    if Annotated_3UTR_file == '':
        print >> sys.stderr, "No annotated 3' UTR file!"
        exit(1)
    if Output_result_file == '':
        print >> sys.stderr, "No result file name!"
        exit(1)
    if sequencing_depth_file == '':
        print >> sys.stderr, "No sequencing depth file!"
        exit(1)

    return Group1_Tophat_aligned_file, Group2_Tophat_aligned_file, output_directory, \
        Annotated_3UTR_file, Output_result_file, Num_least_in_group1_local, \
        Num_least_in_group2_local, Coverage_cutoff_local, Num_threads, sequencing_depth_file


def load_sequencing_depth(depth_file):
    seq_depth_list = []
    for line in open(depth_file, 'r'):
        fields = line.strip('\n').split('\t')
        seq_depth_list.append(int(fields[-1]))
    return np.array(seq_depth_list)


def De_Novo_3UTR_Identification_for_TCGA_Multiple_Samples(argv=None):
    '''shared list with sample id as each element
    '''
    if len(sys.argv) == 1:
        print "Please provide the configure file and specify chr name..."
        exit(1)
    cfg_file = sys.argv[1]
    curr_processing_chr = sys.argv[2]
    print >> sys.stderr, "[%s] Start Analysis ..." % time_now()
    Group1_Tophat_aligned_file, Group2_Tophat_aligned_file, output_directory, \
        Annotated_3UTR_file, Output_result_file, Num_least_in_group1_local, \
        Num_least_in_group2_local, Coverage_cutoff_local, num_threads, sequencing_depth_file \
        = parse_cfgfile(cfg_file)

    num_group_1 = len(Group1_Tophat_aligned_file)
    All_Sample_files = Group1_Tophat_aligned_file[:]
    All_Sample_files.extend(Group2_Tophat_aligned_file)

    global Num_least_in_group1
    global Num_least_in_group2
    global Coverage_cutoff

    if Num_least_in_group1_local != '':
        Num_least_in_group1 = float(Num_least_in_group1_local)
    if Num_least_in_group2_local != '':
        Num_least_in_group2 = float(Num_least_in_group2_local)
    if Coverage_cutoff_local != '':
        Coverage_cutoff = float(Coverage_cutoff_local)

    # Prepare output directory
    output_directory = output_directory.strip('/')+'_'+curr_processing_chr+'/'
    d = os.path.dirname(output_directory)
    if not os.path.exists(d):
        os.makedirs(d)
    temp_dir = d+'/tmp/'
    if not os.path.exists(temp_dir):
        os.makedirs(temp_dir)
    Output_all_prediction_file = output_directory+Output_result_file + '_result_temp.' + curr_processing_chr+'.txt'
    Output_result = open(Output_all_prediction_file, 'w')

    num_samples = len(All_Sample_files)

    print >> sys.stderr, "All tumors Joint Processing %s ..." % curr_processing_chr
    print >> sys.stderr, "[%s] Loading coverage ..." % time_now()

    All_samples_Target_3UTR_coverages, UTR_events_dict = Load_Target_Wig_files_Multiple_threads_shared_dict_sampleid_key(All_Sample_files, Annotated_3UTR_file, num_threads, curr_processing_chr)
    All_samples_sequencing_depths = load_sequencing_depth(sequencing_depth_file)

    print All_samples_sequencing_depths
    All_sample_coverage_weights = All_samples_sequencing_depths/np.mean(All_samples_sequencing_depths)
    # print  All_sample_coverage_weights
    print >> sys.stderr, "[%s] Loading coverage finished ..." % time_now()
    # Write the first line
    first_line = ['Gene', 'fit_value', 'Predicted_Proximal_APA', 'Loci']
    for i in range(num_group_1):
        curr_long_exp = 'A_%s_long_exp' % str(i+1)
        curr_short_exp = 'A_%s_short_exp' % str(i+1)
        curr_ratio = 'A_%s_PDUI' % str(i+1)
        first_line.extend([curr_long_exp, curr_short_exp, curr_ratio])
    for i in range(num_samples - num_group_1):
        curr_long_exp = 'B_%s_long_exp' % str(i+1)
        curr_short_exp = 'B_%s_short_exp' % str(i+1)
        curr_ratio = 'B_%s_PDUI' % str(i+1)
        first_line.extend([curr_long_exp, curr_short_exp, curr_ratio])
    first_line.append('PDUI_Group_diff')
    Output_result.writelines('\t'.join(first_line) + '\n')
    All_events_ids = UTR_events_dict.keys()
    num_threads = 4
    Assigned_events_ids_all_threads = Assign_to_different_processor_balance_events(All_events_ids, num_threads)
    num_real_threads = len(Assigned_events_ids_all_threads)
    Output_each_processor_all = []
    for i in range(num_real_threads):
        curr_temp_output = temp_dir + 'Each_processor_3UTR_Result_%s.txt' % (str(i+1))
        Output_each_processor_all.append(curr_temp_output)

    processs = []
    for i in range(num_real_threads):
        process = multiprocessing.Process(target=Each_Thread_3UTR_estimation_list_version_sample_ids, args=(Assigned_events_ids_all_threads[i], UTR_events_dict, All_sample_coverage_weights, num_samples, num_group_1, Output_each_processor_all[i], All_samples_Target_3UTR_coverages))
        process.start()
        processs.append(process)

    for p in processs:
        p.join()

    # Combine the results
    for i in range(num_real_threads):
        curr_result = Output_each_processor_all[i]
        for line in open(curr_result, 'r'):
            Output_result.writelines(line)
    Output_result.close()
    print >> sys.stderr, "[%s] Finished!" % time_now()


def Each_Thread_3UTR_estimation_list_version_sample_ids(curr_thread_UTR_events_ids, UTR_events_dict, All_sample_coverage_weights, num_samples, num_group_1, Output_result_file, All_samples_coverage_shared_dict):
    Output_result = open(Output_result_file, 'w')
    for curr_3UTR_id in curr_thread_UTR_events_ids:
        curr_3UTR_structure = UTR_events_dict[curr_3UTR_id]
        region_start = curr_3UTR_structure[1]
        region_end = curr_3UTR_structure[2]
        curr_strand = curr_3UTR_structure[-2]
        UTR_pos = curr_3UTR_structure[-1]
        curr_3UTR_all_samples_bp_coverage = []
        for i in range(num_samples):
            curr_sample_curr_3UTR_coverage_wig = All_samples_coverage_shared_dict[curr_3UTR_id, i]
            curr_3UTR_curr_sample_bp_coverage = Convert_wig_into_bp_coverage(curr_sample_curr_3UTR_coverage_wig[0], curr_sample_curr_3UTR_coverage_wig[1], curr_strand)
            curr_3UTR_all_samples_bp_coverage.append(curr_3UTR_curr_sample_bp_coverage)

        select_mean_squared_error, selcted_break_point, UTR_abundances = De_Novo_3UTR_Coverage_estimation_Genome_for_TCGA_multiple_samples(curr_3UTR_all_samples_bp_coverage, region_start, region_end, curr_strand, All_sample_coverage_weights)
        if str(select_mean_squared_error) != "Na":
            num_non_zero = 1
            if num_non_zero > 0:
                All_Long_inclusion_ratios = []
                line_write = [curr_3UTR_id, "%.1f" % select_mean_squared_error, str(selcted_break_point), UTR_pos]
                for i in range(num_samples):
                    if UTR_abundances[0][i] != 'NA':
                        curr_sample_ratio = float(UTR_abundances[0][i])/(float(UTR_abundances[0][i]) + float(UTR_abundances[1][i]))
                        All_Long_inclusion_ratios.append(curr_sample_ratio)
                        # long 3' UTR abundance
                        line_write.append("%.2f" % UTR_abundances[0][i])
                        # short 3' UTR abundance
                        line_write.append("%.2f" % UTR_abundances[1][i])
                        line_write.append("%.2f" % curr_sample_ratio)
                    else:
                        line_write.extend(['NA']*3)

                Output_result.writelines('\t'.join(line_write)+'\n')

    Output_result.close()


def get_version():
    return "0.0.1"


def De_Novo_3UTR_Coverage_estimation_Genome_for_TCGA_multiple_samples(All_Samples_curr_3UTR_coverages, UTR_start, UTR_end, curr_strand, weight_for_second_coverage):
    coverage_threshold = 15
    search_point_start = 150
    search_point_end = int(abs((UTR_end - UTR_start))*0.05)
    num_samples = len(All_Samples_curr_3UTR_coverages)
    # read coverage
    Region_Coverages = []
    Pass_threshold_index = []
    for i in range(num_samples):
        curr_Region_Coverage_raw = All_Samples_curr_3UTR_coverages[i]
        curr_Region_Coverage = curr_Region_Coverage_raw/weight_for_second_coverage[i]
        curr_first_100_coverage = np.mean(curr_Region_Coverage_raw[0:99])
        if curr_first_100_coverage > coverage_threshold:
            Pass_threshold_index.append(i)
            Region_Coverages.append(curr_Region_Coverage)

    least_pass_coverage_num = num_samples*least_pass_coverage_percentage
    if len(Pass_threshold_index) > least_pass_coverage_num and UTR_end - UTR_start >= 150:
        if curr_strand == "+":
            search_region = range(UTR_start+search_point_start, UTR_end-search_point_end+1)
        else:
            search_region = range(UTR_end - search_point_start, UTR_start+search_point_end-1, -1)

        search_region_start = search_point_start
        search_region_end = UTR_end - UTR_start - search_point_end
        Mean_squared_error_list = []
        Estimated_3UTR_abundance_list = []
        for curr_point in range(search_region_start, search_region_end+1):
            curr_search_point = curr_point
            All_samples_result = [[], [], []]
            for curr_sample_region_coverage in Region_Coverages:
                Mean_Squared_error, Long_UTR_abun, Short_UTR_abun = Estimation_abundance(curr_sample_region_coverage, curr_search_point)
                All_samples_result[0].append(Mean_Squared_error)
                All_samples_result[1].append(Long_UTR_abun)
                All_samples_result[2].append(Short_UTR_abun)

            Mean_Squared_error = np.mean(np.array(All_samples_result[0]))
            Mean_squared_error_list.append(Mean_Squared_error)
            Estimated_3UTR_abundance_list.append([All_samples_result[1], All_samples_result[2]])

        if len(Mean_squared_error_list) > 1:
            min_ele_index = Mean_squared_error_list.index(min(Mean_squared_error_list))
            select_mean_squared_error = Mean_squared_error_list[min_ele_index]
            UTR_abundances = Estimated_3UTR_abundance_list[min_ele_index]
            selcted_break_point = search_region[min_ele_index]
            UTR_abundances = [['NA']*num_samples, ['NA']*num_samples]
            UTR_abundances_passed = Estimated_3UTR_abundance_list[min_ele_index]
            for k in range(len(Pass_threshold_index)):
                UTR_abundances[0][Pass_threshold_index[k]] = UTR_abundances_passed[0][k]
                UTR_abundances[1][Pass_threshold_index[k]] = UTR_abundances_passed[1][k]
        else:
            select_mean_squared_error = 'Na'
            UTR_abundances = 'Na'
            selcted_break_point = 'Na'
    else:
        select_mean_squared_error = 'Na'
        UTR_abundances = 'Na'
        selcted_break_point = 'Na'
    return select_mean_squared_error, selcted_break_point, UTR_abundances


def Estimation_abundance(Region_Coverage, break_point):
    Long_UTR_abun = np.mean(Region_Coverage[break_point:])
    Short_UTR_abun = np.mean(Region_Coverage[0:break_point] - Long_UTR_abun)
    if Short_UTR_abun < 0:
        Short_UTR_abun = 0
    Coverage_diff = Region_Coverage[0:break_point] - Long_UTR_abun - Short_UTR_abun
    Coverage_diff = np.append(Coverage_diff, Region_Coverage[break_point:] - Long_UTR_abun)
    Mean_Squared_error = np.mean(Coverage_diff**2)
    return Mean_Squared_error, Long_UTR_abun, Short_UTR_abun


def Load_Target_Wig_files(All_Wig_files, UTR_Annotation_file):
    UTR_events_dict = {}
    All_Samples_Total_depth = []
    for line in open(UTR_Annotation_file, 'r'):
        fields = line.strip('\n').split('\t')
        curr_chr = fields[0]
        region_start = fields[1]
        region_end = fields[2]
        curr_strand = fields[-1]
        UTR_pos = "%s:%s-%s" % (curr_chr, region_start, region_end)
        end_shift = int(round(abs(int(region_start) - int(region_end)) * 0.2))
        if curr_strand == '+':
            region_end = str(int(region_end) - end_shift)
        else:
            region_start = str(int(region_start) + end_shift)
        region_start = int(region_start) + 1
        region_end = int(region_end) - 1
        if region_start + 50 < region_end:
            UTR_events_dict[fields[3]] = [fields[0], region_start, region_end, fields[-1], UTR_pos]

    # Load coverage for all samples
    All_samples_extracted_3UTR_coverage_dict = {}
    for curr_wig_file in All_Wig_files:
        curr_sample_All_chroms_coverage_dict = {}
        num_line = 0
        cur_sample_total_depth = 0
        for line in open(curr_wig_file, 'r'):
            if '#' not in line:
                fields = line.strip('\n').split('\t')
                chrom_name = fields[0]
                region_start = int(fields[1])
                region_end = int(fields[2])
                cur_sample_total_depth += int(fields[-1]) * (region_end - region_start)
                if chrom_name not in curr_sample_All_chroms_coverage_dict:
                    curr_sample_All_chroms_coverage_dict[chrom_name] = [[0], [0]]
                if region_start > curr_sample_All_chroms_coverage_dict[chrom_name][0][-1]:
                    curr_sample_All_chroms_coverage_dict[chrom_name][0].append(region_start)
                    curr_sample_All_chroms_coverage_dict[chrom_name][1].append(0)
                curr_sample_All_chroms_coverage_dict[chrom_name][0].append(region_end)
                curr_sample_All_chroms_coverage_dict[chrom_name][1].append(int(fields[-1]))
            num_line += 1
        curr_sample_All_chroms_coverage_dict[chrom_name][1].append(0)
        All_Samples_Total_depth.append(cur_sample_total_depth)
        for curr_3UTR_event_id in UTR_events_dict:
            curr_3UTR_structure = UTR_events_dict[curr_3UTR_event_id]
            curr_chr = curr_3UTR_structure[0]
            if curr_chr in curr_sample_All_chroms_coverage_dict:
                curr_chr_coverage = curr_sample_All_chroms_coverage_dict[curr_chr]
                region_start = curr_3UTR_structure[1]
                region_end = curr_3UTR_structure[2]
                left_region_index = bisect(curr_chr_coverage[0], region_start)
                right_region_index = bisect(curr_chr_coverage[0], region_end)
                extracted_coverage = curr_chr_coverage[1][left_region_index:right_region_index+1]
                extracted_3UTR_region = curr_chr_coverage[0][left_region_index:right_region_index]
                extracted_3UTR_region.insert(0, region_start)
                extracted_3UTR_region.append(region_end)
                if curr_3UTR_event_id not in All_samples_extracted_3UTR_coverage_dict:
                    All_samples_extracted_3UTR_coverage_dict[curr_3UTR_event_id] = []
                All_samples_extracted_3UTR_coverage_dict[curr_3UTR_event_id].append([extracted_coverage, extracted_3UTR_region])
    return All_samples_extracted_3UTR_coverage_dict, np.array(All_Samples_Total_depth), UTR_events_dict


def Assign_to_different_processor_balance(Total_number, num_processors):
    Assigne_results = []
    num_each_processor = Total_number/num_processors
    if num_each_processor == 0:
        for i in range(Total_number):
            Assigne_results.append([i])
    else:
        remain = Total_number - num_processors * num_each_processor

        for i in range(remain):
            Assigne_results.append(range((i)*(num_each_processor+1), (i+1)*(num_each_processor+1)))

        for i in range(num_processors-remain):
            Assigne_results.append(range(i*num_each_processor+remain*(num_each_processor+1), (i+1)*num_each_processor+remain*(num_each_processor+1)))

    return Assigne_results


def Assign_to_different_processor_balance_events(All_events_ids, num_processors):
    Assigne_results = []
    Total_number = len(All_events_ids)
    num_each_processor = Total_number/num_processors
    if num_each_processor == 0:
        for i in range(Total_number):
            Assigne_results.append([i])

    else:
        remain = Total_number - num_processors * num_each_processor
        for i in range(remain):
            Assigne_results.append(range((i)*(num_each_processor+1), (i+1)*(num_each_processor+1)))

        for i in range(num_processors-remain):
            Assigne_results.append(range(i*num_each_processor+remain*(num_each_processor+1), (i+1)*num_each_processor+remain*(num_each_processor+1)))

    # print Assigne_results
    Assigned_events = []
    print '#assigned events:'
    for curr_processor_inds in Assigne_results:
        curr_processor_events = []
        print len(curr_processor_inds)
        for curr_ele in curr_processor_inds:
            curr_processor_events.append(All_events_ids[curr_ele])
        Assigned_events.append(curr_processor_events)
    return Assigned_events


def Load_Target_Wig_files_Multiple_threads_shared_dict_sampleid_key(All_Wig_files, UTR_Annotation_file, num_threads, curr_processing_chr):
    num_samples = len(All_Wig_files)
    UTR_events_dict = {}
    for line in open(UTR_Annotation_file, 'r'):
        fields = line.strip('\n').split('\t')
        curr_chr = fields[0]
        if curr_chr == curr_processing_chr:
            region_start = fields[1]
            region_end = fields[2]
            curr_strand = fields[-1]
            UTR_pos = "%s:%s-%s" % (curr_chr, region_start, region_end)
            end_shift = int(round(abs(int(region_start) - int(region_end)) * 0.2))
            if curr_strand == '+':
                region_end = str(int(region_end) - end_shift)
            else:
                region_start = str(int(region_start) + end_shift)
            region_start = int(region_start) + 1
            region_end = int(region_end) - 1
            if region_start + 50 < region_end:
                UTR_events_dict[fields[3]] = [fields[0], region_start, region_end, fields[-1], UTR_pos]

    Assigned_index = Assign_to_different_processor_balance(num_samples, num_threads)

    manager = multiprocessing.Manager()
    All_samples_extracted_3UTR_coverage_dict = manager.dict()

    processs = []
    Final_assigned_threads_num = len(Assigned_index)
    for i in range(Final_assigned_threads_num):
        process = multiprocessing.Process(target=load_wig_funct_shared_dict_sampleid_key, args=(All_Wig_files, Assigned_index[i], UTR_events_dict, curr_processing_chr, All_samples_extracted_3UTR_coverage_dict))
        process.start()
        processs.append(process)

    for p in processs:
        p.join()

    return All_samples_extracted_3UTR_coverage_dict, UTR_events_dict


def load_wig_funct_shared_dict_sampleid_key(All_wig_files, assigned_indexes, UTR_events_dict, curr_processing_chr, All_samples_extracted_3UTR_coverage_dict):
    for i in assigned_indexes:
        curr_wig_file = All_wig_files[i]
        print >> sys.stderr, curr_wig_file
        curr_sample_All_chroms_coverage_dict = {}
        with open(curr_wig_file, 'r') as fin:
            for line in fin:
                if line[0] != '#' and line[0] != 't':
                    fields = line.strip('\n').split('\t')
                    chrom_name = fields[0]
                    if chrom_name == curr_processing_chr:
                        region_start = int(fields[1])
                        region_end = int(fields[2])
                        if chrom_name not in curr_sample_All_chroms_coverage_dict:
                            curr_sample_All_chroms_coverage_dict[chrom_name] = [[0], [0]]
                        if region_start > curr_sample_All_chroms_coverage_dict[chrom_name][0][-1]:
                            curr_sample_All_chroms_coverage_dict[chrom_name][0].append(region_start)
                            curr_sample_All_chroms_coverage_dict[chrom_name][1].append(0)
                        curr_sample_All_chroms_coverage_dict[chrom_name][0].append(region_end)
                        curr_sample_All_chroms_coverage_dict[chrom_name][1].append(int(float(fields[-1])))
                    else:
                        if len(curr_sample_All_chroms_coverage_dict) > 0:
                            break
            fin.close()
        if curr_processing_chr not in curr_sample_All_chroms_coverage_dict:
            print >> sys.stderr, 'no wig: '+curr_wig_file
        else:
            curr_sample_All_chroms_coverage_dict[curr_processing_chr][1].append(0)
        curr_sample_coverage_dict = {}
        for curr_3UTR_event_id in UTR_events_dict:
            curr_3UTR_structure = UTR_events_dict[curr_3UTR_event_id]
            curr_chr_local = curr_3UTR_structure[0]
            if curr_chr_local in curr_sample_All_chroms_coverage_dict:
                curr_chr_coverage = curr_sample_All_chroms_coverage_dict[curr_chr_local]
                region_start = curr_3UTR_structure[1]
                region_end = curr_3UTR_structure[2]
                left_region_index = bisect(curr_chr_coverage[0], region_start)
                right_region_index = bisect(curr_chr_coverage[0], region_end)
                extracted_coverage = curr_chr_coverage[1][left_region_index:right_region_index+1]
                extracted_3UTR_region = curr_chr_coverage[0][left_region_index:right_region_index]
                extracted_3UTR_region.insert(0, region_start)
                extracted_3UTR_region.append(region_end)
                curr_event_info = [extracted_coverage, extracted_3UTR_region]
                All_samples_extracted_3UTR_coverage_dict[curr_3UTR_event_id, i] = curr_event_info

Num_least_in_group1 = 1
Num_least_in_group2 = 1
Coverage_cutoff = 30
least_pass_coverage_percentage = 0.3
global All_samples_Target_3UTR_coverages2

if __name__ == '__main__':

    De_Novo_3UTR_Identification_for_TCGA_Multiple_Samples(sys.argv)
