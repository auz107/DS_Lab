import os, sys, csv
sys.path.append('../../')
from tools.userError import userError
from imp import load_source

def load_data_from_tabDelimited_file(inputFile,delType,dataType):
    """
    Imports the data from tab or comma delimited files into a python dictionary 

    INPUTS:
    ------
       inputFile: A string containing the name and path to the input file. The input file should be
                  such that the first row contains the naees of the columns and the first column
                  contains the names of the rows. Example:
     	                col1	col2	col3
                  row1	1.5	1.4	2.5
                  row2	0.1	5.8	4
                  row3	3.5	2.8	5.7

         delType: A string containing the type of delimitation. Eligible choices are 'tab' and 'comma' 
        dataType: A string containing the type of data. Eligible choices are 'float' and 'string'

    OUTPUTS:
    -------
    data: Python dicitonary containing the data

    Ali R. Zomorrodi, Segre Lab @ BU, Last updated 06-19-2015
    """    

    if type(inputFile) is not str:
        raise userError('**Error! inputFile should be a string')
    elif os.path.isfile(inputFile) == False:
        raise userError("**Error! The file '"+inputFile,"' does not exists.")

    if type(delType) is not str:
        raise userError('**Error! delType should be a string')
    elif delType not in ['tab','comma']:
        raise userError("**Error! Invalid delType (the eligible choices are '\t' and ',')")
    elif delType.lower() == 'tab':
        delType = '\t'
    elif delType.lower() == 'comma':
        delType = ','

    if type(dataType) is not str:
        raise userError('**Error! dataType should be a string')
    elif dataType.lower() not in ['float','stirng']:
        raise userError("**Error! Invalid delType (the eligible choices are 'float' and 'string')")

    with open(inputFile,'rb') as file:
        rawData = list(csv.reader(file, delimiter = delType))

    if dataType.lower() == 'float':
        data = dict([((row[0],col),float(row[rawData[0].index(col)])) for row in rawData[1:] for col in rawData[0][1:]])
    elif dataType.lower() == 'string':
        data = dict([((row[0],col),row[rawData[0].index(col)]) for row in rawData[1:] for col in rawData[0][1:]])
    else:
        raise userError('**Error! Invalid dataType even after initial check. Check the code for possible errors')

    return data


def load_data_from_python_file(file_name, var_names):
    """
    Loads data from a file

    INPUTS:
    -------
    file_nmae: 
    Name and paths of the file

    var_names: 
    A list of striings containing the names of the variable in the data filename that must be loaded 

    OUTPUTS:
    --------
    loaded_vars:
    A dictionary where the keys are variables names (those stored in var_names) and values are
    the corresponding loaded variables

    Ali R. Zomorrodi - Segre Lab @ BU
    Last updated: 06/10/2016
    """
    if not isinstance(var_names, list):
        raise TypeError('var_names must be a list')
    elif len([vn for vn in var_names if not isinstance(vn,str)]) > 0:
        raise TypeError('var_names must be a list of string, but non-string objects of type {} were observed in the list'.format(list(set([type(vn) for vn in var_names if not isinstance(vn,str)]))))

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
    
    # Loaded variables
    loaded_vars = dict([(var_name, None) for var_name in var_names])
    for var_name in var_names:
        exec("loaded_vars['" + var_name + "'] = dataFile." + var_name)

    return loaded_vars


