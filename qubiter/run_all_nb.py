import nbformat
from nbconvert.preprocessors import ExecutePreprocessor
import os

'''This script tries to run all jupyter notebooks in the jupyter_notebooks 
folder. The notebooks are executed but not saved (i.e., overwritten) '''

dir_name = 'jupyter_notebooks'
for fname in os.listdir(dir_name):
    if fname[-6:] == '.ipynb':  # and fname[-12:] == 'native.ipynb':
        print("------------", fname)
        try:
            # open() fails when reading markdown Chinese characters
            # and also some types of quotation marks
            nb = nbformat.read(dir_name + '/' + fname,
                               as_version=4)
            ep = ExecutePreprocessor()
            ep.preprocess(nb, {'metadata': {'path': dir_name + "/"}})
        except:
            print('error in ', fname)
        # #this raise will stop execution on first error
        #     raise
        # finally:
        #     nbformat.write(nb, "ERROR_" + fname)
