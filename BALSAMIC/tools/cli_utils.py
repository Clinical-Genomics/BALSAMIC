def add_doc(docstring):
    """
    A decorator for adding docstring. Taken shamelessly from stackexchange. 
    """
    def document(func):
        func.__doc__ = docstring
        return func

    return document
