# Generate Vorticity and Streamfunction Plots from VTK Output of ParMooN


## Installation

### Requirements
 the requirements can be installed from the requirements.txt file using 

 ```
 python -m pip install -r requirements.txt
 ```

 ## Execution

 To run the code, use the following command

 ```
python vorticity.py <vtk file name>
```
 ## Assumptions on Code

 * the code assumes that the u and v velocity components are stored in the vtk file as U_Mean1 and U_Mean2 respectively
