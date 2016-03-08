class userError(Exception):
    """ A calss for returning user defined errors """
    def __init__(self, error = ''):
        if not isinstance(error, str):
            error = '\n**ERROR! ' + str(error) + '\n'
        self.error = '\n**ERROR! ' + error + '\n'

    def __str__(self):
        return self.error


