import os, sys, csv
sys.path.append('../../')
from tools.userError import userError

def importData(inputFile,delType,dataType):
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

