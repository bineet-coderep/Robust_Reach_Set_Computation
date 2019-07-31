# Install Dependencies

Once you have `Python` installed, please install `numpy`.



## Installing Gurobi

* Please obtain appropriate Gurobi License from [here](http://www.gurobi.com/downloads/licenses/license-center). Please note that if you are using Academic License, you **should be in your University network** (VPN should work fine too) while installing the license. Please refer to this [link](https://www.gurobi.com/documentation/8.1/quickstart_windows/academic_validation.html) for details. After the license is installed properly, Gurobi can be used from home network.
* Install Gurobi. Please note that we will need Gurobi Python Interface. On-line documentation on installation can be found [here](http://www.gurobi.com/documentation/).
* Gurobi Python Interface can also be installed through [Anaconda](https://www.anaconda.com/). Details on installing Gurobi Python Interface through `conda` can be found [here](https://www.gurobi.com/documentation/8.1/quickstart_mac/installing_the_anaconda_py.html#section:Anaconda).



## Testing the Environment

_Please note that this step is optional_.

The above installation can be checked by running [this]https://github.com/bineet-coderep/Robust_Reach_Set_Computation/blob/master/LinearSystemsWithFaults/EnvironmentTester/test-install.py) python script in `LinearSystemWithFaults/EnvironmentTester/test-install.py`

​	`python test-install.py`

If all the installations were correctly done, then this script should run without any error and display `Environment is Ready` at the end.
