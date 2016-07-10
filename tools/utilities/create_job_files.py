from __future__ import division
from datetime import timedelta  # To convert elapsed time to hh:mm:ss format
import re

def create_intervals(total_cases_num,interval_size, one_case_req_time = None, stdout_msgs = True):
    """
    Creates the intervals dividing a large interval into smaller ones

    INPUTS:
    -------
    total_cases_num: 
    The total number of cases to consider 

    interval_size: 
    Desired size of the iteration intervals

    one_case_req_time:
    The required time in seconds to run only case of the problem. 
    If this parameter is provided it is used to estimate how long 
    each thread (interval) would take

    OUTPUTS:
    --------
    intervals:
    Intervals sum of which giving total_cases_num

    Ali R. Zomorrodi - Segre Lab @ BU
    Last updated: 06/03/2016 
    """
    intervals = []

    # Start positions in each interval

    start_positions = range(1,total_cases_num + 1,interval_size)
    for k in start_positions:
            if k + (interval_size - 1) <= total_cases_num:
                intervals.append((k,k + (interval_size - 1)))
            else:
                intervals.append((k,total_cases_num))

    if one_case_req_time != None:
        thread_time_s = interval_size*one_case_req_time 
        thread_time_h = thread_time_s/3600
        thread_time_d = thread_time_h/24

    if stdout_msgs:
        print '\n'
        print '**Total # of intervals = ',len(intervals),'\n'

        for interval in intervals:
           print interval

        if one_case_req_time != None:
            print '\nEach thread would take {} (days, )hours:minutes:seconds to complete\n'.format(str(timedelta(seconds = thread_time_s)))

    return intervals

def create_job_files(total_cases_num, interval_size, job_filename_base, joboutput_filename_base, results_filename_base, results_filename_input_format, start_pos_input_format, end_pos_input_format, main_code_dir, commands, max_walltime = 24):
    """
    Creates the job files given:

    INPUTS:
    -------
    total_cases_num: 
    Total number of combinations 

    interval_size: 
    The desired of iteration intervals

    outfile_base_name: 
    The base name of the output files

    python_commad: 
    The task that should be performed ('analyze_games' or 'analyze_cost')

    job_filename_base: 
    Base name of the output job file

    joboutput_filename_base: 
    The name of the file storing the dynamic output of the job. It is essentially the same as job_filename_base
    with the file path deleted (in order to create .out files in the same directory as the job file)

    results_filename_base: 
    Base name of the results file for the main python script

    results_filename_input_format:
    This is a string showing how the results_filename is provided in the input of the function.
    For example, if a funciton is in the form func(x,y, results_file), 
    results_filename_input_format = 'results_file = None'. So, always provide the keyword 
    for results_filename foloowed by one space followed by None. Then this script replaces None
    in with the created results_filename (from results_filename_base) in any command in commands
    that it appears 

    start_pos_input_format & end_post_input_format:
    The same as results_filename_input_format for the start and end position of the loops

    main_code_dir:
    A string containing the directory where the main code run by the job file is located

    commands:
    A list of strings containing the commands that msut be written to files. If the command are in the form
    'pyton -c "..."', it's best to use triple-quoted string format after python -c
    """
    intervals = create_intervals(total_cases_num = total_cases_num,interval_size = interval_size)

    for interval in intervals:
        job_filename = job_filename_base + '_' + str(interval[0]) + '_' + str(interval[1]) + '.sh'
        joboutput_filename = joboutput_filename_base + '_' + str(interval[0]) + '_' + str(interval[1]) + '.out'
        print 'creating file ',job_filename,'...'
        with open(job_filename,'w') as outfile:
            outfile.write('#!/bin/bash -l\n')
            outfile.write('#$-l h_rt=' + str(int(max_walltime))+ ':00:00\n\n')
            outfile.write('cd ' + main_code_dir + '\n\n')

            # This is to merge the .e[jobid] and .o[jobid] files
            outfile.write('#$ -j y\n')

            # Set the output file
            outfile.write('\n#$ -o ' + joboutput_filename + '\n')

            outfile.write('\nsource activate /projectnb/bioinfor/SEGRE/alizom\n\n')

            # python -c "import time;print '\n**Job started at ',time.strftime('%c'),'\n'" > job_dynamic_yeast_stsrt_pos_end_pos.out 2>&1
            outfile.write("\npython -c \"import time;print '\\n**Job started at ',time.strftime('%c'),'\\n'\"\n")

            start_pos = interval[0]
            end_pos = interval[1]
            results_filename = results_filename_base + "_" + str(start_pos) + "_" + str(end_pos) + '.py'

            for command in commands:
                start_pos_input = re.sub(' = None', ' = ' + str(start_pos), start_pos_input_format)
                command = re.sub(start_pos_input_format, start_pos_input, command)

                end_pos_input = re.sub(' = None', ' = ' + str(end_pos), end_pos_input_format)
                command = re.sub(end_pos_input_format, end_pos_input, command)

                results_filename_input = re.sub(' = None', " = '" + results_filename + "'", results_filename_input_format)
                command = re.sub(results_filename_input_format, results_filename_input, command)

                outfile.write('\n' + command +  '\n')

            # python -c "import time;print '\n**Job ended at ',time.strftime('%c'),'\n'" >> job_dynamic_yeast_start_pos_end_pos.out 2>&1
            outfile.write("\npython -c \"import time;print '\\n**Job ended at ',time.strftime('%c'),'\\n'\"\n")

