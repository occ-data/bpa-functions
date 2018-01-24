USC BloodPAC Data Commons Query 1
December 4, 2017

Purpose:
    •Reproduce a simple query across experiments on the USC dataset 
    •Get some level of QA to make sure the USC submission looks as expected   
    •Use these methods to generate meaningful analysis at a slightly deeper data level

Query:
This query will show the level of depth of analysis for the samples analyzed.  The query will retrieve a list of all patients
 (de-identified patient IDs) and determine which assays were performed on each patient sample.  

Implementation:
    •Run test in the data portal for prototyping the query using GraphQL (data.bloodpac.org)
    •Build out some functions in python in the BloodPAC Virtual Machine
    •Perform query in the VM using python and generate visualization of resulting data 

Jupyter Notebook File:
USC_Patient_EXPs.ipynb


Output:
The plot shows the results of the query with the de-identified Patient ID on the vertical axis and the type of assay performed
 on the horizontal axis. A green line all the way across represents a patient sample that has all 6 assays performed.


================================================================================
================================================================================
 
 
 
USC BloodPAC Data Commons Query 2
January 17, 2018

Purpose:
    •Recreate graphs used in a publication based off submitted data
    •Get some level of QA to make sure the USC submission looks as expected   

Query:
This query will create 3 Swarmplots comparing DNA Concentration by cell classication in respect to 24 or 48 hour times to fractionizations.

Implementation:
    •Run test in the data portal for prototyping the query using GraphQL (data.bloodpac.org)
    •Build out some functions in python in the BloodPAC Virtual Machine
    •Perform query in the VM using python and generate visualization of resulting data 

Jupyter Notebook File:
USC_CTC_Plots.ipynb


Output:
The plot shows the results of the query side by side in respect to cell classication ("CTC", "CTC Low", "CTC Small"). The x-axis is either
24 or 48 hours to fractionization. The y-axis is the measurement in DNA concentration.

