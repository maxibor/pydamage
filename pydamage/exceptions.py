class Error(Exception):
    """Base class for exceptions in this module."""
    pass


class AlignmentFileError(Error):
    """Exception raised for errors in validating the dataset agains the standards.

    Attributes:
        message -- explanation of the error
    """

    def __init__(self, message):
        self.message = message
      