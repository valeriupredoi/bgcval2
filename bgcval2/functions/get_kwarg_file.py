import os

def get_kwarg_file(kwargs, filekey, default = False):
    """
    From the kwarg dict and the filekey, this function 
    find the filepath. 

    It doesn't make any effort to load the file though, just the path 
    as a string
    """
    # load file, area and mask
    filename = kwargs.get(filekey, default)

    if not filename:
        raise FileNotFoundError(f"Function requires an {filekey} kwarg to run calculation.")

    if isinstance(filename, list) and len(filename)>1:
        raise KeyError(f'Several {filekey} files found: {filename}')

    if isinstance(filename, list) and len(filename)==1:
        filename = filename[0]

    if  not os.path.exists(filename):
        raise FileNotFoundError(f"Function {filekey} kwarg file not found: {filename}")

    return filename

