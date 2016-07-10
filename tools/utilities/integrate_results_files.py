from load_data_fromFile import load_data_from_python_file
    
def integrate_results_files(results_filenames, results_varname, output_filename):
    """
    Integrates the results of all scripts (when splitting one job to multiple 
    smaller jobs) inot one sinble file 
 
    INPUTS:
    -------
    results_filenames: 
    A list of strings containing the names of the files containing the results

    results_varname: 
    Name of the variable storing the results

    output_filename: 
    Name of the output file name containing the integration of all results

    Ali R. Zomorrodi - Daniel Segre Lab @ BU
    Last updated: 06/16/2016
    """
    # sum of the entries in all result files
    entries_sum = 0

    results = {}
    # Import the data in the module stored in file_name
    for file_name in results_filenames:
        results_curr = load_data_from_python_file(file_name = file_name, var_names = [results_varname])[results_varname]

        entries_sum += len(results_curr.keys())
        print 'Total # of entries in ', file_name,' = ',len(results_curr.keys())

        results_keys = results.keys()
        for k in results_curr.keys():
            if k in results_keys:
                raise userError(str(k) + ' in ' + file_name + ' already in results_keys\n')
            else:
                results[k] = results_curr[k]

    print '\nThe total # of entries in results = {},  entries_sum = {} '.format(len(results.keys()),entries_sum)

    # Write the results in the specified output file
    with open(output_filename,'w') as f:
        f.write(results_varname + ' = {}\n')
        for k in results.keys():
            f.write(results_varname + '[' + str(k) + '] = ' + str(results[k]) + '\n')
    print '\nResults were integrated and written into ',output_filename,'\n'

def del_results_files(results_filenames):
    """
    Deletes the results files after verification

    INPUTS:
    -------
    results_filenames: A list of strings containing the names of the files containing the results
    """
    for file_name in results_filenames:
        print 'deleteing ',file_name, ' ...'
        os.system('rm -rf ' + file_name)

#----------------------
if __name__ == '__main__':
    pass
