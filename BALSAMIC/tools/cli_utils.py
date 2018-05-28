import os

def add_doc(docstring):
    """
    A decorator for adding docstring. Taken shamelessly from stackexchange. 
    """
    def document(func):
        func.__doc__ = docstring
        return func

    return document

def createDir(path):
    if os.path.isdir(path):
        basepath = os.path.basename(path)
        basepath_number = 0
        if "." in basepath:
            basepath_number = int(basepath.split(".")[1])
        basepath_string = basepath.split(".")[0]
        basepath_number += 1
        path = os.path.join(os.path.dirname(path), ".".join([basepath_string,str(basepath_number)]))
        createDir(path)
    else:
        os.makedirs(os.path.dirname(path), exist_ok=True)
        return True
