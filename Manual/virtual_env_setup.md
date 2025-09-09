# Setting up a virtual environment

### What and why?
One of the strengths of Python is that there are lots of great packages available to leverage when scripting up a project. But Python is not great at managing these packages, which can lead to your Python installation getting cluttered and conflicts arising when packages are updated. Virtual environments help deal with this problem by setting up a stand--alone folder structure with all the packages needed for a specific project installed in an isolated Python environment. Virtual environments are particularly useful for keeping consistency between users when the code is being shared. Conda is used for managing environments and package installations with the Anaconda distribution of Python. 

### Setting up a Conda environment
A file containing the Conda environment settings for Bryan is included with the source code:  *_env_bryan_.yml*. To import the environment:

1.	Copy the file *_env_bryan.yml* to a local folder (e.g. C:\Users\username\Documents)
2.	Install Anaconda distribution (if not already installed)
3.	From the start menu, type and click on Anaconda Prompt.

This should open up a command line interface for Anaconda. In this command line interface undertake the following:

1. Use the *cd* command to navigate to the folder where your local copy of *_env_bryan.yml* is located.
2. Enter the command: conda env create -f _env_bryan.yml
3. Verify that the environment was set up successfully by entering: conda activate bryan29
 
The bryan29 environment is activated from the batch file that is used to run Bryan. An example batch file is shown below.

```dos
:: Insert the name of the main config file
set config_file="C:\path\to\sims_config.json"

:: Insert the file of the python project
set pyfile="C:\path\to\Bryan\Main.py"

:: Activate the Anaconda bryan29 environment
call C:\Users\username\Anaconda3\envs\bryan29\activate.bat

:: Run the model
python %pyfile% %config_file%

pause
```