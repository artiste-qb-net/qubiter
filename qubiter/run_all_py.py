import os
from subprocess import Popen, PIPE, call

'''This script tries to run all class files, in all folders. Class files all 
have a main() at the end. '''


dir_whitelist = [
    "./",
    'adv_applications',
    'device_specific',
    'latex_tools',
    'quantum_CSD_compiler'
    ]
file_blacklist = [
    '__init__.py',
    'run_all_nb.py',
    'run_all_py.py',
    'class_diagram.py'
]
for dir_name in dir_whitelist:
    for fname in os.listdir(dir_name):
        if fname[-3:] == '.py' and fname not in file_blacklist:
            path = dir_name + '/' + fname
            print('--------------------', path)
            pro = Popen(['python', fname], cwd=dir_name,
                        stdout=PIPE, stderr=PIPE)
            stdout, stderr = pro.communicate()
            print(str(stderr))

