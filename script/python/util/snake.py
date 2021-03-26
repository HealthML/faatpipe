
'''
functions and classes that we want to be accessible in python scripts and snakemake
'''

def clean_str(x):
    '''
    function that removes problematic symbols from strings wo we can use them as file-prefixes
    '''
    return x.replace('(','').replace(')','')
