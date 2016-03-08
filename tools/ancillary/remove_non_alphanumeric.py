import re

def remove_non_alphanumeric(input_string, replace_with_underline = False):
    """
    This functions replaces the following non-alphanumeric characters in the input string
    or the input list of string input_string:
    if replace_with_underline is True
        Replaces + with _plus_
        replaces everything else with underline

    if replace_with_underline is False
        Replaces + with _plus_
        Replaces - with _dash_
        Replaces ' with _prime_
        Replaces " with _doublePrime_
        Replaces ( with _LParen_
        Replaces ) with _RParen_
        Replaces , with _comma_
        Replaces all other non-alphanumeric characters with underline

    INPUTS:
    ------
    input_string: A string or a list of strings
    replace_with_underline: If True, dahs, coma, left and right paranetheses are replaced with an underline score
    OUTPUTS:
    a string with all non-alaphabetical and non-numerical characters replaced with an underline

    Ali R. Zomorrodi - Segre Lab @ Boston University
    Last updated: 01-12-2016 
    """

    converted_string = []

    # If the input is a string convert it to a list with one element
    if isinstance(input_string,str):
        input_string = [input_string]

    elif not isinstance(input_string,list):
        raise userError('Invluad input for remove_non_alphanumeric! String or list of strings expected.')

    for s in input_string:
        # Remove space at the end and at the begining 
        s = re.sub(' $','',s)
        s = re.sub('^ ','',s)

        # Replace - with _dash
        if replace_with_underline:
            s = re.sub('\-','_',s)
        else:
            s = re.sub('\-','_dash_',s)

        # Replace + with _plus_
        s = re.sub('\+','_plus_',s)

        # Replace "'" with _prime 
        if replace_with_underline:
            s = re.sub("'",'_',s)
        else:
            s = re.sub("'",'_prime_',s)

        # Replace '"' with _doublePrime 
        if replace_with_underline:
            s = re.sub('"','_',s)
        else:
            s = re.sub('"','_doublePrime_',s)

        # Replace "(" with "_LParen_" 
        if replace_with_underline:
            s = re.sub('\(','_',s)
        else:
            s = re.sub('\(','_LParen_',s)

        # Replace ")" with "_RParen_" 
        if replace_with_underline:
            s = re.sub('\)','_',s)
        else:
            s = re.sub('\)','_RParen_',s)

        # Replace "," with "_comma" 
        if replace_with_underline:
            s = re.sub(',','_',s)
        else:
            s = re.sub(',','_comma_',s)

        # Replace any other non-alphanumeric character with an underline
        s = re.sub('[^\w]','_',s)

        # Replace any two subsequent underlines with one
        while '__' in s:
            s = re.sub('__','_',s)

        # Repmove any underline at the begining or at the end of a name
        s = re.sub('^_','',s)
        s = re.sub('_$','',s)

        converted_string.append(s)

    if len(converted_string) == 1:
        converted_string = converted_string[0]

    return converted_string

