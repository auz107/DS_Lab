class userError(Exception):
    """ 
    A calss for returning user defined errors 

    See the section titled "8.5. User-defined Exceptions" in Python's documentation:
    https://docs.python.org/2/tutorial/errors.html
    """
    def __init__(self, error_msgs = ''):
        self.error = '\n**ERROR! ' + str(error_msgs) + '\n'

    def __str__(self):
        return self.error


