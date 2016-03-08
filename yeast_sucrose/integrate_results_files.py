import os
from imp import load_source
    
def integrate_results_files(results_filenames, results_var_name, output_file_name):
    """
    Integrates the results of all scripts (when splitting one job to multiple smaller jobs) inot one sinble file 
 
    INPUTS:
    -------
    results_filenames: A list of strings containing the names of the files containing the results
      results_var_name: Name of the variable storing the results
      output_file_name: Name of the output file name containing the integration of all results
    """
    # sum of the entries in all result files
    entries_sum = 0

    results = {}
    # Import the data in the module stored in file_name
    for file_name in results_filenames:
        if type(file_name) == str:
            if not os.path.isfile(file_name):
                raise IOError("No such file was found :'" + file_name + "'")
            else:
                # First delete the model dataFile if it already exists. If it is not deleted
                # the new module data is merged with the previous ones
                try:
                    del sys.modules['dataFile']
                except:
                    pass
                load_source('dataFile',file_name)
                import dataFile
        exec('results_curr = dataFile.' + results_var_name)
        
        entries_sum += len(results_curr.keys())
        print 'Total # of entries in ', file_name,' = ',len(results_curr.keys())

        results_keys = results.keys()
        for k in results_curr.keys():
            if k in results_keys:
                raise userError(str(k) + ' in ' + file_name + ' already in results_keys\n')
            else:
                results[k] = results_curr[k]

    print '\nThe total # of entries in results = {},  entries_sum = {} '.format(len(results.keys()),entries_sum)

    # Write the results inot the specified output file
    with open(output_file_name,'w') as f:
        f.write(results_var_name + ' = {}\n')
        for k in results.keys():
            f.write(results_var_name + '[' + str(k) + '] = ' + str(results[k]) + '\n')
    print '\nResults were integrated and written into ',output_file_name,'\n'

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
    print '\n---- Cost results, his, fixed o2 ---'
    results_filenames = ['results/cost_results_his_fixedo2_1_600.py','results/cost_results_his_fixedo2_601_1071.py'] 
    integrate_results_files(results_filenames = results_filenames, results_var_name = 'cooperation_cost', output_file_name = 'results/cost_results_his_fixedo2.py')
    #del_results_files(results_filenames = results_filenames)

    print '\n---- Cost results, his, independent o2 ---'
    results_filenames = ['results/cost_results_his_indepO2_1_600.py','results/cost_results_his_indepO2_601_1071.py'] 
    integrate_results_files(results_filenames = results_filenames, results_var_name = 'cooperation_cost', output_file_name = 'results/cost_results_his_indepO2.py')
    #del_results_files(results_filenames = results_filenames)

    print '\n---- Cost results, his, variable o2 ---'
    results_filenames = ['results/cost_results_his_variableO2_1_600.py','results/cost_results_his_variableO2_601_1071.py'] 
    integrate_results_files(results_filenames = results_filenames, results_var_name = 'cooperation_cost', output_file_name = 'results/cost_results_his_variableO2.py')
    #del_results_files(results_filenames = results_filenames)

    print '\n---- Cost results, atp, fixed o2 ---'
    results_filenames = ['results/cost_results_atp_fixedo2_1_500.py', 'results/cost_results_atp_fixedo2_501_1000.py', 'results/cost_results_atp_fixedo2_1001_1500.py', 'results/cost_results_atp_fixedo2_1501_2000.py', 'results/cost_results_atp_fixedo2_2001_2500.py', 'results/cost_results_atp_fixedo2_2501_3000.py', 'results/cost_results_atp_fixedo2_3001_3500.py', 'results/cost_results_atp_fixedo2_3501_4000.py', 'results/cost_results_atp_fixedo2_4001_4500.py', 'results/cost_results_atp_fixedo2_4501_5000.py', 'results/cost_results_atp_fixedo2_5001_5151.py'] 
    integrate_results_files(results_filenames = results_filenames, results_var_name = 'cooperation_cost', output_file_name = 'results/cost_results_atp_fixedo2.py')
    #del_results_files(results_filenames = results_filenames)

    print '\n---- Cost results, atp, independent o2 ---'
    results_filenames = ['results/cost_results_atp_indepO2_1_500.py', 'results/cost_results_atp_indepO2_501_1000.py', 'results/cost_results_atp_indepO2_1001_1500.py', 'results/cost_results_atp_indepO2_1501_2000.py', 'results/cost_results_atp_indepO2_2001_2500.py', 'results/cost_results_atp_indepO2_2501_3000.py', 'results/cost_results_atp_indepO2_3001_3500.py', 'results/cost_results_atp_indepO2_3501_4000.py', 'results/cost_results_atp_indepO2_4001_4500.py', 'results/cost_results_atp_indepO2_4501_5000.py', 'results/cost_results_atp_indepO2_5001_5151.py'] 
    integrate_results_files(results_filenames = results_filenames, results_var_name = 'cooperation_cost', output_file_name = 'results/cost_results_atp_indepO2.py')
    #del_results_files(results_filenames = results_filenames)

    print '\n---- Cost results, atp, variable o2 ---'
    results_filenames = ['results/cost_results_atp_variableO2_1_500.py', 'results/cost_results_atp_variableO2_501_1000.py', 'results/cost_results_atp_variableO2_1001_1500.py', 'results/cost_results_atp_variableO2_1501_2000.py', 'results/cost_results_atp_variableO2_2001_2500.py', 'results/cost_results_atp_variableO2_2501_3000.py', 'results/cost_results_atp_variableO2_3001_3500.py', 'results/cost_results_atp_variableO2_3501_4000.py', 'results/cost_results_atp_variableO2_4001_4500.py', 'results/cost_results_atp_variableO2_4501_5000.py', 'results/cost_results_atp_variableO2_5001_5151.py'] 
    integrate_results_files(results_filenames = results_filenames, results_var_name = 'cooperation_cost', output_file_name = 'results/cost_results_atp_variableO2.py')
    #del_results_files(results_filenames = results_filenames)

    print '\n---- Game results, his, fixed o2 ---'
    results_filenames = ['results/ss_game_results_his_fixedo2_1_600.py','results/ss_game_results_his_fixedo2_601_1071.py'] 
    integrate_results_files(results_filenames = results_filenames, results_var_name = 'results', output_file_name = 'results/ss_game_results_his_fixedo2.py')
    #del_results_files(results_filenames = results_filenames)

    print '\n---- Game results, his, independent o2 ---'
    results_filenames = ['results/ss_game_results_his_indepO2_1_600.py','results/ss_game_results_his_indepO2_601_1071.py'] 
    integrate_results_files(results_filenames = results_filenames, results_var_name = 'results', output_file_name = 'results/ss_game_results_his_indepO2.py')
    #del_results_files(results_filenames = results_filenames)

    print '\n---- Games results, his, variable o2 ---'
    results_filenames = ['results/ss_game_results_his_variableO2_1_600.py','results/ss_game_results_his_variableO2_601_1071.py'] 
    integrate_results_files(results_filenames = results_filenames, results_var_name = 'results', output_file_name = 'results/ss_game_results_his_variableO2.py')
    #del_results_files(results_filenames = results_filenames)

    print '\n---- Games results, atp, fixed o2 ---'
    results_filenames = ['results/ss_game_results_atp_fixedo2_1_500.py', 'results/ss_game_results_atp_fixedo2_501_1000.py', 'results/ss_game_results_atp_fixedo2_1001_1500.py', 'results/ss_game_results_atp_fixedo2_1501_2000.py', 'results/ss_game_results_atp_fixedo2_2001_2500.py', 'results/ss_game_results_atp_fixedo2_2501_3000.py', 'results/ss_game_results_atp_fixedo2_3001_3500.py', 'results/ss_game_results_atp_fixedo2_3501_4000.py', 'results/ss_game_results_atp_fixedo2_4001_4500.py', 'results/ss_game_results_atp_fixedo2_4501_5000.py', 'results/ss_game_results_atp_fixedo2_5001_5151.py'] 
    integrate_results_files(results_filenames = results_filenames, results_var_name = 'results', output_file_name = 'results/ss_game_results_atp_fixedo2.py')
    #del_results_files(results_filenames = results_filenames)

    print '\n---- Games results, atp, independent o2 ---'
    results_filenames = ['results/ss_game_results_atp_indepO2_1_500.py', 'results/ss_game_results_atp_indepO2_501_1000.py', 'results/ss_game_results_atp_indepO2_1001_1500.py', 'results/ss_game_results_atp_indepO2_1501_2000.py', 'results/ss_game_results_atp_indepO2_2001_2500.py', 'results/ss_game_results_atp_indepO2_2501_3000.py', 'results/ss_game_results_atp_indepO2_3001_3500.py', 'results/ss_game_results_atp_indepO2_3501_4000.py', 'results/ss_game_results_atp_indepO2_4001_4500.py', 'results/ss_game_results_atp_indepO2_4501_5000.py', 'results/ss_game_results_atp_indepO2_5001_5151.py'] 
    integrate_results_files(results_filenames = results_filenames, results_var_name = 'results', output_file_name = 'results/ss_game_results_atp_indepO2.py')
    #del_results_files(results_filenames = results_filenames)

    print '\n---- Game results, atp, variable o2 ---'
    results_filenames = ['results/ss_game_results_atp_variableO2_1_500.py', 'results/ss_game_results_atp_variableO2_501_1000.py', 'results/ss_game_results_atp_variableO2_1001_1500.py', 'results/ss_game_results_atp_variableO2_1501_2000.py', 'results/ss_game_results_atp_variableO2_2001_2500.py', 'results/ss_game_results_atp_variableO2_2501_3000.py', 'results/ss_game_results_atp_variableO2_3001_3500.py', 'results/ss_game_results_atp_variableO2_3501_4000.py', 'results/ss_game_results_atp_variableO2_4001_4500.py', 'results/ss_game_results_atp_variableO2_4501_5000.py', 'results/ss_game_results_atp_variableO2_5001_5151.py'] 
    integrate_results_files(results_filenames = results_filenames, results_var_name = 'results', output_file_name = 'results/ss_game_results_atp_variableO2.py')
    #del_results_files(results_filenames = results_filenames)


