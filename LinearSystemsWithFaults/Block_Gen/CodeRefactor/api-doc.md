# API Documentation - Uncertain Linear System Verifier



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
* `isFeasible():` Returns `True` if the **Hard Coded** unsafe condition is satisfied by the reached set given by A^k starting from the initial set as given in `initialSet`. 



## Class Benchmarks

_This class is used to represent various benchmarks_

Each class representing a Benchmark should have following fields:

 * `A`: A `numpy` square array representing A of the given dynamics of the form Ax+B.
 * `B`: A `numpy` 2-D array representing B of the given dynamics of the form Ax+B.
 * `mode`: It can take values from the set `{'+','.'}`. `'+'` represents the given dynamics is discrete and `'.'` represents the given dynamics is continuous. 
 * `h`: If the given dynamics is continuous, then `h` is the step size used for discretizing it. 



## Class BenchmarkVerify

### Methods

* `verifyTemplate(Matrix A,Matrix B,int n1,int n2,int h,char mode,list f,list v,int T,list initset,int uS,float uC,String strName)`:  Displays the _'Final Verification Report'_  (Explained in `how-to-use.md`) of the dynamics with the given parameters. Following are explanation of the parameters:
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

* `readFromFile(String fname)`: Displays the _'Final Verification Report'_ of the dynamics as given in the input XML file `fname`. Format of the input file is explained in `how-to-use.md`
* **Other methods are also available for Inbuilt Benchmark Verification**



  

