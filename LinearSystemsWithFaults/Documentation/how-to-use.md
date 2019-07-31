# How to Use

This document is on how to use the main verification script `uncertain-system-verifier.py` ([link](https://github.com/bineet-coderep/Robust_Reach_Set_Computation/blob/master/LinearSystemsWithFaults/VerificationEngine/uncertain-system-verifier.py)) in `LinearSystemsWithFaults/VerificationEngine/`

This script can take input by following methods:

* XML File (Two different formats for constant and time-varying perturbations)
* Inbuilt Benchmarks
* Define your dynamics as `sub classes` in `Class Benchmarks`. 

This script can be used to generate the following output:

* Verification Report (format explained later) of the given input dynamics with the given faults (Constant or Time-Varying).
* Verification Report (format explained later) of the given input dynamics without any fault.
* Possible places/blocks which satisfies the conditions mentioned in the paper.

The Verification Report contains the following information:

*  Safety status up-to the given time `T` _i.e._ Safe or Unsafe.
* Total Time taken by the Verification Engine.
* If the system is unsafe, it can also display a counter example. (if this option is chosen)

Please note that customizations (or extensions) can be done using the available APIs (please refer to the API Documentation at `LinearSystemsWithFaults/Documentation/api-doc.md` ([here](https://github.com/bineet-coderep/Robust_Reach_Set_Computation/blob/refactoring/LinearSystemsWithFaults/Documentation/api-doc.md)))





## Input



* ### XML File (Constant Perturbation)

The template of the input XML format can be found in `LinearSystemsWithFaults/VerificationEngine/InputFiles/test-input-template.xml` ([here](https://github.com/bineet-coderep/Robust_Reach_Set_Computation/blob/refactoring/LinearSystemsWithFaults/VerificationEngine/InputFiles/test-input-template.xml))

The whole input should be enclosed with the tag `<input>`

A sample XML input file can be found in `LinearSystemsWithFaults/VerificationEngine/InputFiles/test-input2.xml` ([here](https://github.com/bineet-coderep/Robust_Reach_Set_Computation/blob/refactoring/LinearSystemsWithFaults/VerificationEngine/InputFiles/test-input2.xml))

Following are the description of other tags:

#### Tag `<A>`

This should contain the matrix A of the dynamics of the form Ax+B

Each line within this tag represents a row; and at each line the values should be separated by `,`  

#### Tag `<B>`

This should contain the matrix B of the dynamics of the form Ax+B

Each line within this tag represents a row; and at each line the values should be separated by `,`  

#### Tag `<h>` 

If the given dynamics is continuous, then `h` is the step size used for discretizing it. 

This should be single value within this tag

#### Tag `<mode>` 

Denotes the dynamics is discrete or continuous.

This should be single character within this tag. It can take values from the set `{'+','.'}`. `'+'` represents the given dynamics is discrete and `'.'` represents the given dynamics is continuous. 

#### Tag `<faults>`

 Denotes the set of blocks with uncertainties

Each line within this tag represents a block. And at each line four values are separated by `,` representing the tuple `(startrow,length,startcol,breadth)`

#### Tag `<variance>`

An ordered list of tuple which denotes the range, the uncertain variable at that index can vary.

Each line with this tag represent the range of uncertainty for the block at the same line number. And at each line, two values are separated using `,` denoting `(min,max)` range.

#### Tag `<initialset>`

An ordered list of tuple of the form `(int min, int max)` representing the range of value for the state variable at that index.

Each line with this tag represent the range of uncertainty for the block at the same line number. And at each line, two values are separated using `,` denoting `(min,max)` range.

#### Tag `<name>`

Name of the model.

The value at this tag should be string.

#### Tag `<time>`

Time up-to which verification is desired.

The value at this tag should be an integer.

#### Tag `<unsafestate>`

The state number for which the unsafe condition is to be checked.

The value at this tag should be an integer.

#### Tag `<unsafeparam>`

The unsafe state `uS` should be less than equals `uC`.

The value at this tag should be an integer.



* ### XML File (Time-Varying Perturbation)

The template of the input XML format can be found in `LinearSystemsWithFaults/VerificationEngine/InputFiles/test-input-templateTV.xml` ([here](https://github.com/bineet-coderep/Robust_Reach_Set_Computation/blob/refactoring/LinearSystemsWithFaults/VerificationEngine/InputFiles/test-input-templateTV.xml))

The whole input should be enclosed with the tag `<input>`

A sample XML input file can be found in `LinearSystemsWithFaults/VerificationEngine/InputFiles/test-inputTV2.xml` ([here](https://github.com/bineet-coderep/Robust_Reach_Set_Computation/blob/refactoring/LinearSystemsWithFaults/VerificationEngine/InputFiles/test-inputTV2.xml))

Following are the description of other tags:

#### Tag `<A>`

This should contain the matrix A of the dynamics of the form Ax+B

Each line within this tag represents a row; and at each line the values should be separated by `,`  

#### Tag `<B>`

This should contain the matrix B of the dynamics of the form Ax+B

Each line within this tag represents a row; and at each line the values should be separated by `,`  

#### Tag `<h>` 

If the given dynamics is continuous, then `h` is the step size used for discretizing it. 

This should be single value within this tag

#### Tag `<mode>` 

Denotes the dynamics is discrete or continuous.

This should be single character within this tag. It can take values from the set `{'+','.'}`. `'+'` represents the given dynamics is discrete and `'.'` represents the given dynamics is continuous. 

#### Tag `<faults>`

 Denotes the set of blocks with uncertainties

Each line within this tag represents a block. And at each line four values are separated by `,` representing the tuple `(startrow,length,startcol,breadth)`

#### Tag `<variance>`

An ordered list of list of tuple which denotes the range, the uncertain variable at that index can vary.

Each line with this tag represent the range of time-varying uncertainty for the block at the same line number. And at each line, the range of uncertainty for each time step is separated by ';'. And for each step two values are separated using `,` denoting `(min,max)` range.

Eg: At line number l: `min-t0,max-t0; min-t1,max-t1; min-t2,max-t2`. 

The range of perturbation corresponding to the uncertain variable `l`, at time step 1 is `(min-t1,max-t1)`. 

#### Tag `<initialset>`

An ordered list of tuple of the form `(int min, int max)` representing the range of value for the state variable at that index.

Each line with this tag represent the range of uncertainty for the block at the same line number. And at each line, two values are separated using `,` denoting `(min,max)` range.

#### Tag `<name>`

Name of the model.

The value at this tag should be string.

#### Tag `<time>`

Time up-to which verification is desired.

The value at this tag should be an integer.

#### Tag `<unsafestate>`

The state number for which the unsafe condition is to be checked.

The value at this tag should be an integer.

#### Tag `<unsafeparam>`

The unsafe state `uS` should be less than equals `uC`.

The value at this tag should be an integer.

#### Tag `<step>`

The number of steps after which the perturbation range varies as per the input given in tag `<variance>`.



* ### Inbuilt Benchmarks

A list of Benchmarks are already defined inside the `Class Benchmarks` (Please see `api-doc.md` for details).

These set of predefined benchmarks are used inside `Class Verify`, `Class VerifyFaultFree` and `Class ConsolidatedVerificationReport`. Inside these classes inbuilt functions are available (`def <benchmark-name>`), where you can modify the parameters (like faults, initial set, time at etc) as desired.



* ### Define your own Dynamics as a Sub-Class

Define a sub-class inside the `Class Benchmarks` (Please see `api-doc.md` for details)

Then define a function (just like the predefined ones of the form `def <benchmark-name`) inside `Class Verify`, `Class VerifyFaultFree` or/and `Class ConsolidatedVerificationReport` (Also in `class Search` if you want to find possible faulty blocks). Inside this function you can define your parameters like faults, initial set, time and etc.





## Running the Script for Verification

* ### XML as Input File (Constant Perturbation)

  To get Verification Report of the given system (as XML file) with and without fault together, run the following command:

  ​	`python uncertain-system-verifier.py <fname>`

  ​				or

  ​	`ConsolidatedVerificationReport.readFromFile("<file-name>")` as API function call

  

  To get Verification Report of the given system (as XML file) with fault, call the following function:

  ​	`Verify.readFromFile("<file-name>")`

  

  To get Verification Report of the given system (as XML file) without fault, call the following function:

  ​	`VerifyFaultFree.readFromFile("<file-name>")`

  

* ### XML as Input File (Time-Varying Perturbation)

  To get Verification Report of the given system (as XML file) with Time-Varying fault, call the following function:

  ​	`Verify.readFromFileTimeVarying("<file-name>")`

  

  To get Verification Report of the given system (as XML file) with and without fault together, call the following function:

  ​	`ConsolidatedVerificationReport.readFromFileTimeVarying("<file-name>")`

  

* ### Inbuilt Benchmarks

  To get Verification Report of desired benchmark (Inbuilt) with constant fault and without fault together, call the following function:

  ​	`ConsolidatedVerificationReport.<benchmark-name>()`
  For example: `ConsolidatedVerificationReport.flightEnvelope()`

  To get Verification Report of desired benchmark (Inbuilt) with time-varying fault and without fault together, call the following function:

  ​	`ConsolidatedVerificationReport.<benchmark-name>TimeVarying()`
  For example: `ConsolidatedVerificationReport.flightEnvelopeTimeVarying()`

  

  To get Verification Report of desired benchmark (Inbuilt) with constant fault, call the following function:

  ​	`Verify.<benchmark-name>()`
  For example: `Verify.flightEnvelope()`

  To get Verification Report of desired benchmark (Inbuilt) with time-varying fault, call the following function:

  ​	`Verify.<benchmark-name>TimeVarying()`
  For example: `Verify.flightEnvelopeTimeVarying()`

  

  To get Verification Report of desired benchmark (Inbuilt) without fault, call the following function:

  ​	`VerifyFaultFree.<benchmark-name>()`
  For example: `VerifyFaultFree.flightEnvelope()`

  

* ### User Defined Dynamics

  To get Verification Report of a user defined system with and without fault together, call the following function:

  ​		`ConsolidatedVerificationReport.<system-name>()` or 
  

​		`ConsolidatedVerificationReport.<system-name>Timevarying()` 
  For example: `ConsolidatedVerificationReport.test()` or `ConsolidatedVerificationReport.testTimeVarying()` 

  

To get Verification Report of defined system with fault, call the following function:

​	`Verify.<system-name>()`or `Verify.<system-name>TimeVarying()` 
  For example: `Verify.test()` or `Verify.testTimeVarying()`



  To get Verification Report of defined system without fault, call the following function:

  ​	`VerifyFaultFree.<benchmark-name>()`
  For example: `VerifyFaultFree.test()`






## Running the Script for Finding Possible Faults

* ### XML as Input File

  To get the list of blocks of various sizes in the given dynamics (as XML), that satisfies the two conditions mentioned in the paper, call the following function:

  ​	`Search.readFromFile("<file-name>")`

* ### Inbuilt Benchmarks

  To get the list of blocks of various sizes in the desired Benchmark dynamics (Inbuilt), that satisfies the two conditions mentioned in the paper, call the following function:

  ​	`Search.<benchmark-name>()`
  For example:
  ​	`Search.flightEnvelope()`

* ### User Defined Dynamics

  To get the list of blocks of various sizes in the defined dynamics, that satisfies the two conditions mentioned in the paper, call the following function:

  ​	`Search.<system-name>()`
  For example:
  ​	`Search.madeUp()`
