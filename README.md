# bpa-functions
This repository contains tools for the BloodPAC community to answer and ask questions about the BloodPAC data commons. 

# Organization

The demo notebook and bpa_analysis_functions_v2.py are provided to help the community get started understanding and using the BloodPAC data commons.   

## Python: bpa_analysis_functions

The python functions library contains commands to accomplish basic tasks like:

* Adding API keys to your profile
* Getting files from a bucket
* Querying the BloodPAC API

It also contains a variety of custom functions to accomplish tasks like:

* Query and summarize fields
* Create plots from summarized fields
* Create lists of files
* etc.

To use add a line `import bpa_analysis_functions_v2.py as bpa`. 

## Community Notebooks

To further build the BloodPAC community, members are contributing notebooks around their own analysis.   To add your work, please clone the repo, add your notebook, then submit a pull request.  

Community notebooks can also be used to answer questions being asked in the BloodPAC Q&A forum.  

https://groups.google.com/d/forum/bloodpac-q-and-a

NOTE:  since this is a public repository, be sure that any notebooks you post obfuscate any sensitive data results.   If you're running a Jupyter Notebook that has sensitive information, please be sure to submit a pull request with the cells not run. 
