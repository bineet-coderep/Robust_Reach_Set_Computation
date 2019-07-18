# API Documentation - Uncertain Linear System Verifier

The API python file can be found in `/LinearSystemsWithFaults/VerificationEngine/RobustReachSetAPI.py` ([link](https://github.com/bineet-coderep/Robust_Reach_Set_Computation/blob/refactoring/LinearSystemsWithFaults//VerificationEngine/RobustReachSetAPI.py))

## Class Structure

### Constructor

* ```Structure(list l,int n)```: Where ```list l``` is a list of tuple of the form ``` (startrow,length,startcol,breadth)``` representing the places with 1 in the support. And ```int n``` is the dimension of the support matrix.

### Methods

* ```visualizeStr()```: Displays the support matrix.



## Class BlockStructure

_Inherits class ``Structure``_

### Constructor

* ```BlockStructure(ind, int n)```: Where ```ind``` is of the form ```(startrow,length,startcol,breadth)``` representing the block containing 1. And ```int n``` is the dimension of the support matrix. 

### Methods

* ```getRowStart()```: Returns the starting row of the Block (of 1s).
* ```getRowLength()```: Returns the number of rows in the Block.
* ```getColStart()```:  Returns the starting column of the Block (of 1s).
* ```getColBreadth```: Returns the number of columns in the Block.
* ```mutualKernelStrCheck(BlockStructure b)```: Returns ```True``` if the product of this structure with b is 0. ```False``` otherwise.
* ```getFaultyCells()```: Returns the `list` of cells contained in the block.



## Class StructLME

### Constructor

* `StructLME(list l)`: Where `list l` is a list of `Structure`. It represents the support matrices of all the uncertain variables (N_i) and the constant support too (N_0).

### Methods

* `isCondSAT()`: Returns `True`  if the `list` of `Structure` supports the two sufficient conditions mentioned in the paper.
* `checkStrN0()`: Returns `True` if the support of N_0 satisfies the conditions mentioned in the paper.
* `visualizeStructLME()`: Displays the list of support matrices.
* `findStrN0ForOne(BlockStructure blk)`: _(Static Method)_ Returns a super support for N_0 such that the given block support `blk` satisfies the conditions mentioned in the paper.
* `intersectionOfBlocks(Structure X, Structure Y, int n)`: Returns a `Structure` which is intersection of `X` and `Y` of size `n`.



## Class Dynamics

### Constructor

* `Dynamics(A, int n, list v)`: `A` is the matrix representing the linear dynamics, of size `n*n` and `v` is the list of faulty blocks represented as ```(startrow,length,startcol,breadth)```.

### Methods

* `displayA()`: Displays the matrix `A` representing the linear dynamics.
* `findAllFaultyCells()`: Returns a list of cells which are uncertain/faulty in the given dynamics represented by `A` .
* `findConstStr()`: Returns the Support Matrix (`Structure` of N_0) of the given Dynamics (Represented by the matrix `A`).
* `findConstMat()`: Return the N_0 matrix corresponding to the LME of dynamics `A`.
* `findVarBlkStr()`: Returns the list block support matrices representing the faulty variables of the dynamics (as given in `v`).
* `findVarBlkMats():`  Returns the list of N_i (for i>=1) matrices corresponding to the LME of dynamics `A`.



## Class DynamicsRepOp

### Constructor

* `DynamicsRepOp(N0, list N, int s)`: The Linear Matrix Expression of the dynamics is given as parameter. N0 is the matrix corresponding to N_0 and the `list N` corresponds to the other N_i matrices (for i>=1). And `int s` is the dimension of the matrices N_i.

### Methods

* `displayA()`: Displays the Linear Matrix Expression.
* `power(int k)`: Calculates the A^k as given in the paper; where A is the matrix corresponding to the LME. Returns a `list` representing the LME.
* `power(int k)`: Calculates the A^k as given in the paper; where A is the matrix corresponding to the LME. Returns a `list` representing the LME.
* `powerTimeVarying(int k, int n)`: Calculates the A^k as given in the paper for time-varying perturbation; where A is the matrix corresponding to the LME and `n` is the number of steps after which the perturbation range changes . Returns a `list` representing the LME.



## Class InputFormation

### Constructor

* `InputFormation(Matrix A, int n, list vars)`: `A` is the matrix of size `n*n` representing the dynamics. `list vars` is the list of blocks with uncertainties, represented using the tuple `(startrow,length,startcol,breadth)`. 

### Methods

* `formMatA()`: Changes the value of all the faulty cells in `A` to 1.
* `getAllFaultyCells()`: Returns the list of cells which correspond to some faulty variable.



## Class CalcReachableSet

### Constructor

* `CalcReachableSet(DynamicsReOp A, int k, list initialSet, list fault, list perturbation, int unsafe-state, float unsafe-param)`: `list initialSet` is an ordered list of tuple of the form `(int min, int max)` representing the range of value for the state variable at that index. `list fault` is the list of blocks with uncertainties, represented using the tuple `(startrow,length,startcol,breadth)`. `list perturbation` is  a ordered list of tuple of the form `(int min, int max)` representing the range of values it can take, corresponding to the faulty variable (`list fault`) at the same index. The last two parameters correspond to unsafe condition check; _i.e._ the state number `unsafe-state` is less than equals the value `unsafe-param`

### Methods

* `getAllFaultyCells()`: Returns the list of cells which correspond to some faulty variable.
* `indexofBlk(i,j)`: Returns the index of the block support (represented in the field `faults`) corresponding to the cell `(i,j)`.
* `isFeasible():` Returns `True` if the unsafe condition (the value of the `unsafe-state` th state variable is less than equals `unsafe-param`) is satisfied by the reached set given by A^k starting from the initial set as given in `initialSet`. 



## Class CalcReachableSetTimeVarying

### Constructor

- `CalcReachableSetTimeVarying(DynamicsReOp A, int k, list initialSet, list fault, list perturbation, int unsafe-state, float unsafe-param, int t)`: `list initialSet` is an ordered list of tuple of the form `(int min, int max)` representing the range of value for the state variable at that index. `list fault` is the list of blocks with uncertainties, represented using the tuple `(startrow,length,startcol,breadth)`. `list perturbation` is  an ordered list of list of tuple of the form `(int min, int max)` representing the range of values it can take at each time step, corresponding to the faulty variable (`list fault`) at the same index. The last two parameters correspond to unsafe condition check; _i.e._ the state number `unsafe-state` is less than equals the value `unsafe-param`. `t` is number of step after which the perturbation varies.

### Methods

- `getAllFaultyCells()`: Returns the list of cells which correspond to some faulty variable.
- `indexofBlk(i,j)`: Returns the index of the block support (represented in the field `faults`) corresponding to the cell `(i,j)`.
- `isFeasible():` Returns `True` if the unsafe condition (the value of the `unsafe-state` th state variable is less than equals `unsafe-param`) is satisfied by the reached set given by A^k starting from the initial set as given in `initialSet`. 
- `isFeasible:` Returns `True` if the unsafe condition (the value of the `unsafe-state` th state variable is less than equals `unsafe-param`) is satisfied by the reached set given by A^k starting from the initial set as given in `initialSet` with the time-varying faults as given by `perturbation`.



## Class CalcReachableSetFaultFree

### Constructor

- `CalcReachableSetFaulFree(DynamicsReOp A, int n, int k, list initialSet, int unsafe-state, float unsafe-param)`: `int n` is the dimension of `A`, `list initialSet` is an ordered list of tuple of the form `(int min, int max)` representing the range of value for the state variable at that index. The last two parameters correspond to unsafe condition check; _i.e._ the state number `unsafe-state` is less than equals the value `unsafe-param`

### Methods

- `isFeasible():` Returns `True` if the unsafe condition (the value of the `unsafe-state` th state variable is less than equals `unsafe-param`) is satisfied by the reached set given by A^k starting from the initial set as given in `initialSet`. 



## Class Benchmarks

_This class is used to represent various benchmarks_

Each class representing a Benchmark should have following fields:

 * `A`: A `numpy` square array representing A of the given dynamics of the form Ax+B.
 * `B`: A `numpy` 2-D array representing B of the given dynamics of the form Ax+B.
 * `mode`: It can take values from the set `{'+','.'}`. `'+'` represents the given dynamics is discrete and `'.'` represents the given dynamics is continuous. 
 * `h`: If the given dynamics is continuous, then `h` is the step size used for discretizing it. 



## Class Verify

### Methods

* `verifyTemplate(Matrix A,Matrix B,int h,char mode,list f,list v,int T,list initset,int uS,float uC,String strName)`:  Displays the _'Final Verification Report'_  (Explained in `how-to-use.md`) of the dynamics with the given parameters. Following are explanation of the parameters:
  * `Matrix A`: A `numpy` square array representing A of the given dynamics of the form Ax+B.
  * `Matrix B`: A `numpy` 2-D array representing B of the given dynamics of the form Ax+B.
  * `char mode`: It can take values from the set `{'+','.'}`. `'+'` represents the given dynamics is discrete and `'.'` represents the given dynamics is continuous. 
  * `int h`: If the given dynamics is continuous, then `h` is the step size used for discretizing it. 
  * `list f`: List of blocks with uncertainties, represented with the tuple `(startrow,length,startcol,breadth)` 
  * `list v`: An ordered list of integers which denotes the percentage, the uncertain variable at that index can vary.
  * `int T`: Time up to which verification of the system is desired.
  * `list initset`: An ordered list of tuple of the form `(int min, int max)` representing the range of value for the state variable at that index.
  * `int uS`: The state number for which the unsafe condition is to be checked.
  * `float uC`: The unsafe state `uS` should be less than equals `uC`.
  * `String strName`: Name of the Benchmark.
* `verifyTemplateNonVerbose(Matrix A,Matrix B,int h,char mode,list f,list v,int T,list initset,int uS,float uC)`:  Returns a dictionary with the _'Final Verification Report'_  (Explained in `how-to-use.md`) of the dynamics with the given parameters. Following are explanation of the parameters:
  - `Matrix A`: A `numpy` square array representing A of the given dynamics of the form Ax+B.
  - `Matrix B`: A `numpy` 2-D array representing B of the given dynamics of the form Ax+B.
  - `char mode`: It can take values from the set `{'+','.'}`. `'+'` represents the given dynamics is discrete and `'.'` represents the given dynamics is continuous. 
  - `int h`: If the given dynamics is continuous, then `h` is the step size used for discretizing it. 
  - `list f`: List of blocks with uncertainties, represented with the tuple `(startrow,length,startcol,breadth)` 
  - `list v`: An ordered list of tuple which denotes the range, the uncertain variable at that index can vary.
  - `int T`: Time up to which verification of the system is desired.
  - `list initset`: An ordered list of tuple of the form `(int min, int max)` representing the range of value for the state variable at that index.
  - `int uS`: The state number for which the unsafe condition is to be checked.
  - `float uC`: The unsafe state `uS` should be less than equals `uC`.
* `verifyTemplateTimeVarying(Matrix A,Matrix B,int h,char mode,list f,list v,int T,list initset,int uS,float uC,String strName, int s)`:  Displays the _'Final Verification Report'_  (Explained in `how-to-use.md`) of the dynamics with the given parameters for time-varying perturbation. Following are explanation of the parameters:
  - `Matrix A`: A `numpy` square array representing A of the given dynamics of the form Ax+B.
  - `Matrix B`: A `numpy` 2-D array representing B of the given dynamics of the form Ax+B.
  - `char mode`: It can take values from the set `{'+','.'}`. `'+'` represents the given dynamics is discrete and `'.'` represents the given dynamics is continuous. 
  - `int h`: If the given dynamics is continuous, then `h` is the step size used for discretizing it. 
  - `list f`: List of blocks with uncertainties, represented with the tuple `(startrow,length,startcol,breadth)` 
  - `list v`: An ordered list of list of tuple which denotes the range of perturbations at every time-step for each uncertain variable. 
  - `int T`: Time up to which verification of the system is desired.
  - `list initset`: An ordered list of tuple of the form `(int min, int max)` representing the range of value for the state variable at that index.
  - `int uS`: The state number for which the unsafe condition is to be checked.
  - `float uC`: The unsafe state `uS` should be less than equals `uC`.
  - `String strName`: Name of the Benchmark.
  - `int s`: The perturbation will vary after `s` steps according to `list v`.
* `verifyTemplateTimeVaryingNonVerbose(Matrix A,Matrix B,int h,char mode,list f,list v,int T,list initset,int uS,float uC,String strName, int s)`:  Returns a dictionary with the _'Final Verification Report'_  (Explained in `how-to-use.md`) of the dynamics with the given parameters for time-varying perturbations. Following are explanation of the parameters:
  - `Matrix A`: A `numpy` square array representing A of the given dynamics of the form Ax+B.
  - `Matrix B`: A `numpy` 2-D array representing B of the given dynamics of the form Ax+B.
  - `char mode`: It can take values from the set `{'+','.'}`. `'+'` represents the given dynamics is discrete and `'.'` represents the given dynamics is continuous. 
  - `int h`: If the given dynamics is continuous, then `h` is the step size used for discretizing it. 
  - `list f`: List of blocks with uncertainties, represented with the tuple `(startrow,length,startcol,breadth)` 
  - `list v`: An ordered list of list of tuple which denotes the range of perturbations at every time-step for each uncertain variable. 
  - `int T`: Time up to which verification of the system is desired.
  - `list initset`: An ordered list of tuple of the form `(int min, int max)` representing the range of value for the state variable at that index.
  - `int uS`: The state number for which the unsafe condition is to be checked.
  - `float uC`: The unsafe state `uS` should be less than equals `uC`.
  - `String strName`: Name of the Benchmark.
  - `int s`: The perturbation will vary after `s` steps according to `list v`.
* `readFromFile(String fname)`: Displays the _'Final Verification Report'_ of the dynamics as given in the input XML file `fname`. Format of the input file is explained in `how-to-use.md`.
* `readFromFileTimeVarying(String fname)`: Displays the _'Final Verification Report'_ of the dynamics as given in the input XML file `fname` for time-varying perturbation. Format of the input file is explained in `how-to-use.md`
* **Other methods are also available inbuilt, for various Benchmark Verification**. Inbuilt methods are available for both constant and time-varying noise as well. Method names with <benchmark> are to verify the benchmark with a constant perturbation whereas Method names with <benchmark> followed by `TimeVarying` are to verify with time-varying noise.



## Class VerifyFaultFree

### Methods

- `verifyTemplateFaultFree(Matrix A,Matrix B,int h,char mode,int T,list initset,int uS,float uC,String strName)`:  Displays the _'Final Verification Report'_  without any uncertainty (Explained in `how-to-use.md`) of the dynamics with the given parameters. Following are explanation of the parameters:
  - `Matrix A`: A `numpy` square array representing A of the given dynamics of the form Ax+B.
  - `Matrix B`: A `numpy` 2-D array representing B of the given dynamics of the form Ax+B.
  - `char mode`: It can take values from the set `{'+','.'}`. `'+'` represents the given dynamics is discrete and `'.'` represents the given dynamics is continuous. 
  - `int h`: If the given dynamics is continuous, then `h` is the step size used for discretizing it. 
  - `int T`: Time up to which verification of the system is desired.
  - `list initset`: An ordered list of tuple of the form `(int min, int max)` representing the range of value for the state variable at that index.
  - `int uS`: The state number for which the unsafe condition is to be checked.
  - `float uC`: The unsafe state `uS` should be less than equals `uC`.
  - `String strName`: Name of the Benchmark.
- `verifyTemplateFaultFreeNonVerbose(Matrix A,Matrix B,int h,char mode,int T,list initset,int uS,float uC)`:  Returns a dictionary with the _'Final Verification Report'_  without any uncertainty (Explained in `how-to-use.md`) of the dynamics with the given parameters. Following are explanation of the parameters:
  - `Matrix A`: A `numpy` square array representing A of the given dynamics of the form Ax+B.
  - `Matrix B`: A `numpy` 2-D array representing B of the given dynamics of the form Ax+B.
  - `char mode`: It can take values from the set `{'+','.'}`. `'+'` represents the given dynamics is discrete and `'.'` represents the given dynamics is continuous. 
  - `int h`: If the given dynamics is continuous, then `h` is the step size used for discretizing it. 
  - `int T`: Time up to which verification of the system is desired.
  - `list initset`: An ordered list of tuple of the form `(int min, int max)` representing the range of value for the state variable at that index.
  - `int uS`: The state number for which the unsafe condition is to be checked.
  - `float uC`: The unsafe state `uS` should be less than equals `uC`.
- `readFromFile(String fname)`: Displays the _'Final Verification Report'_ without any uncertainty of the dynamics as given in the input XML file `fname`. Format of the input file is explained in `how-to-use.md`
- **Other methods are also available Inbuilt for various Benchmark Verification**



## Class Finder

### Constructor

* `Finder(Matrix A, Matrix B, float h, char mode)`: `Matrix A` and `Matrix B`  defines a dynamics of the form Ax+B, where `mode` is either `+` (Discrete) or `.` (Continuous). And `float h` is the discretization step size.

### Methods

* `enumAllBlocks()`: Returns a list enumerating all possible blocks (of all sizes) that satisfies the sufficient conditions mentioned in the paper. 
* `viewStructLME()`: Displays the structure of LMEs and if it satisfies the conditions mentioned in the paper.
* `enumBlock(m,k)`: Returns a list enumerating all possible blocks of size `m*k` that satisfies the sufficient conditions mentioned in the paper.



## Class Search

### Methods

* `searchTemplate(Matrix A, Matrix B, float h, char mode)`:  Displays the _'Blocks Satisfying the Conditions mentioned in the Paper'_  (Explained in `how-to-use.md`) of the dynamics with the given parameters. Following are explanation of the parameters:
  - `Matrix A`: A `numpy` square array representing A of the given dynamics of the form Ax+B.
  - `Matrix B`: A `numpy` 2-D array representing B of the given dynamics of the form Ax+B.
  - `char mode`: It can take values from the set `{'+','.'}`. `'+'` represents the given dynamics is discrete and `'.'` represents the given dynamics is continuous. 
  - `int h`: If the given dynamics is continuous, then `h` is the step size used for discretizing it. 
* `readFromFile(String fname)`: Displays the _'Blocks Satisfying the Conditions mentioned in the Paper'_  of the dynamics as given in the input XML file `fname`. Format of the input file is explained in `how-to-use.md`
* **Other methods are also available inbuilt for Finding Blocks in various Benchmarks** 



## Class ConsolidatedVerificationReport

### Methods

* `conVerTemplate(Matrix A,Matrix B,int h,char mode,list f,list v,int T,list initset,int uS,float uC,String strName)`: Displays the _'Final Verification Report'_  with constant uncertainty and without uncertainty (Explained in `how-to-use.md`) of the dynamics with the given parameters. Following are explanation of the parameters:
  - `Matrix A`: A `numpy` square array representing A of the given dynamics of the form Ax+B.
  - `Matrix B`: A `numpy` 2-D array representing B of the given dynamics of the form Ax+B.
  - `char mode`: It can take values from the set `{'+','.'}`. `'+'` represents the given dynamics is discrete and `'.'` represents the given dynamics is continuous. 
  - `int h`: If the given dynamics is continuous, then `h` is the step size used for discretizing it. 
  - `int T`: Time up to which verification of the system is desired.
  - `list initset`: An ordered list of tuple of the form `(int min, int max)` representing the range of value for the state variable at that index.
  - `int uS`: The state number for which the unsafe condition is to be checked.
  - `float uC`: The unsafe state `uS` should be less than equals `uC`.
  - `String strName`: Name of the Benchmark.
* `conVerTemplateTimeVarying(Matrix A,Matrix B,int h,char mode,list f,list v,int T,list initset,int uS,float uC,String strName, int step)`: Displays the _'Final Verification Report'_  with Time-Varying uncertainty and without uncertainty (Explained in `how-to-use.md`) of the dynamics with the given parameters. Following are explanation of the parameters:
  - `Matrix A`: A `numpy` square array representing A of the given dynamics of the form Ax+B.
  - `Matrix B`: A `numpy` 2-D array representing B of the given dynamics of the form Ax+B.
  - `char mode`: It can take values from the set `{'+','.'}`. `'+'` represents the given dynamics is discrete and `'.'` represents the given dynamics is continuous. 
  - `int h`: If the given dynamics is continuous, then `h` is the step size used for discretizing it. 
  - `int T`: Time up to which verification of the system is desired.
  - `list initset`: An ordered list of tuple of the form `(int min, int max)` representing the range of value for the state variable at that index.
  - `int uS`: The state number for which the unsafe condition is to be checked.
  - `float uC`: The unsafe state `uS` should be less than equals `uC`.
  - `String strName`: Name of the Benchmark.
  - `int Step`: Number of steps after which the range of perturbation changes.
* `readFromFile(String fname)`: Displays the _'Final Verification Report'_  with constant uncertainty and without any uncertainty of the dynamics as given in the input XML file `fname`. Format of the input file is explained in `how-to-use.md`
* `readFromFileTimevarying(String fname)`: Displays the _'Final Verification Report'_  with Time-varying uncertainty and without any uncertainty of the dynamics as given in the input XML file `fname`. Format of the input file is explained in `how-to-use.md`
* **Other methods are also available Inbuilt for various Benchmark Verification (with both Time-varying and constant perturbations)**





  

