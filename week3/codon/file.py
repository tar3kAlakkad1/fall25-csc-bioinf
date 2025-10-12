__all__ = ["InvalidFileError", "TreeError"]


class TreeError(Static[Exception]):
    """
    An exception that occurs in context of tree topology.
    """

    error_message: str

    def __init__(self, error_message: str):
        super().__init__("TreeError", "Exception occured in context of tree topology.")
        self.error_message = error_message


class InvalidFileError(Static[Exception]):
    """
    Indicates that the file is not suitable for the requested action,
    either because the file does not contain the required data or
    because the file is malformed.
    """

    error_message: str

    def __init__(self, error_message: str):
        super().__init__(
            "InvalidFileError",
            "the file is not suitable for the requested action, either because the file does not contain the required data or because the file is malformed. ",
        )
        self.error_message = error_message
