def parse_notes(note_string, note_delimiter = ':'):
    """
    This is function parse_legacy_sbml_notes from sbml.py in COBRApy 
    """
    note_dict = {}
    start_tag = '<p>'
    end_tag = '</p>'
    if '<html:p>' in note_string:
        start_tag = '<html:p>'
        end_tag = '</html:p>'
    while start_tag in note_string and end_tag in note_string:
        note_start = note_string.index(start_tag)
        note_end = note_string.index(end_tag)
        the_note = note_string[(note_start + len(start_tag)):note_end].lstrip(' ').rstrip(' ')
        if note_delimiter in the_note:
            note_delimiter_index = the_note.index(note_delimiter)
            note_field = the_note[:note_delimiter_index].lstrip(' ').rstrip(' ').replace('_',' ').upper()
            note_value = the_note[(note_delimiter_index+1):].lstrip(' ').rstrip(' ')
            if note_field in note_dict:
                note_dict[note_field ].append(note_value)
            else:
                note_dict[note_field] = [note_value]
        note_string = note_string[(note_end+len(end_tag)): ]

    if 'CHARGE' in note_dict and note_dict['CHARGE'][0].lower() in ['none', 'na', 'nan']:
        note_dict.pop('CHARGE') #Remove non-numeric charges
        
    return(note_dict)

