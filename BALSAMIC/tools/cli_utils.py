import os


def add_doc(docstring):
    """
    A decorator for adding docstring. Taken shamelessly from stackexchange. 
    """

    def document(func):
        func.__doc__ = docstring
        return func

    return document


def createDir(path, interm_path=[]):
    if os.path.isdir(os.path.abspath(path)):
        basepath = os.path.basename(os.path.abspath(path))
        basepath_number = 0
        if "." in basepath:
            basepath_number = int(basepath.split(".")[1])
        basepath_string = basepath.split(".")[0]
        basepath_number += 1
        path = os.path.join(
            os.path.dirname(os.path.abspath(path)),
            ".".join([basepath_string, str(basepath_number)]))
        interm_path.append(path)
        createDir(path, interm_path)
        return interm_path[-1]
    else:
        os.makedirs(os.path.abspath(path), exist_ok=True)
        return interm_path[-1]
