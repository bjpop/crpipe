'''
Various utility functions that don't have a sensible home elsewhere
'''

import os

def safe_make_dir(path):
    '''Make a directory if it does not already exist'''
    if not os.path.exists(path):
        os.makedirs(path)
