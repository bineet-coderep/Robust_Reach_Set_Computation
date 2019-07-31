# Robust Reachable Set Computation

Based on the paper: 

_Robust Reachable Set: Accounting for Uncertainties in Linear Dynamical Systems_ by Bineet Ghosh and Parasara Sridhar Duggirala. (_EMSOFT 2019, To appear_)




## Dependencies

This Verification Engine has following dependencies:

* Python Numpy.
* [Gurobi](https://www.gurobi.com/) Python Interface.

The installation of the above dependencies can be checked by running [this](https://github.com/bineet-coderep/Robust_Reach_Set_Computation/blob/master/LinearSystemsWithFaults/EnvironmentTester/test-install.py) python script, as follows:

​		`python test-install.py`

If all the packages were installed correctly, this script should run without any error and display `Environment is Ready!!`

Details on installation can be found in `/LinearSystemsWithFaults/Docuementation/install.md` ([link](https://github.com/bineet-coderep/Robust_Reach_Set_Computation/blob/master/LinearSystemsWithFaults/Documentation/install.md))



## How to Use

The main verification script is `/LinearSystemsWithFaults/VerificationEngine/uncertain-system-verifier.py`.

Details on how to use it to verify linear dynamical systems with uncertainties or verify inbuilt benchmarks (also modifying them) is documented in `/LinearSystemsWithFaults/Docuementation/how-to-use.md` ([link](https://github.com/bineet-coderep/Robust_Reach_Set_Computation/blob/master/LinearSystemsWithFaults/Documentation/how-to-use.md))

For customizing the behavior of the input and output, or using more features please refer to the API Documentation in `/LinearSystemsWithFaults/Docuementation/api-doc.md` ([link](https://github.com/bineet-coderep/Robust_Reach_Set_Computation/blob/master/LinearSystemsWithFaults/Documentation/api-doc.md))



## API Documentation

The API documentation can be found in `/LinearSystemsWithFaults/Docuementation/api-doc.md` ([link](https://github.com/bineet-coderep/Robust_Reach_Set_Computation/blob/master/LinearSystemsWithFaults/Documentation/api-doc.md))

Please send an email to Bineet Ghosh at ghosh.bineet22@gmail.com if you have any questions.



## Bug Report

Please send an email to Bineet Ghosh at ghosh.bineet22@gmail.com to report any bug.
