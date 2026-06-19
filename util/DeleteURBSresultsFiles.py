"""
This script removes all *.p and *.il files in a folder and subfolders
"""

import glob
import os

path_name = 'C:/PythonProjects/TFD_2024/03_Design/runs/E005_repeat/US_E005_repeat/urbs/model/**/'
extensions = ['*.p', '*.il']

for extension in extensions:
    full_path = f'{path_name}{extension}'
    files_in_dir = glob.iglob(full_path, recursive = True)   
    for _file in files_in_dir:
        print(_file) 
        os.remove(_file)