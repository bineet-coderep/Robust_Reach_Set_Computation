Title         : Reachable Set Computation for Uncertain Linear Systems Implementation

[TITLE]

# `Class Structure`
_Represents Structure of a matrix as defined in definition 3.1_

## Constructor
`Structure(l,n)`: A list containing a tuple of indices (startrow, length, startcol, breadth). And n is the size of the matrix.

## Fields

  * `omegas`: A list containing a tuple of indices (startrow, length, startcol, breadth).

  * `size`: Dimension of the matrix.

  * `matrix`: Creates an actual matrix representing this structure. 1 represents the set R, and 0 otherwise.

## Functions
  * `visualizeStr()`: Displays the matrix

# `Class BlockStructure`
_Inherits the class Structure. Represents block structured matrices as per definition 3.3_

## Constructor
`BlockStructure(l,n)`: Where l is tuple of (startrow, length, startcol, breadth) and size is n

## Functions
  * `getRowStart()`: Returns he starting row of the block structure

  * `getRowLength()`: Returns the number of rows spanned by the block structure

  * `getColStart()`: Returns the starting column of the block structure

  * `getColBreadth()`: Returns the number of columns spanned by the block structure

  * `intersection(x,lx,y,ly)`: A static method. Returns true iff the intervals [[x,lx]] and [[y,ly]] intersects.

  * `mutualKerStrCheck(mat)`: Returns True iff mat is Mutual Kernel Structure (Definition 3.4) of the calling object.

# `Class RelationGraph`
_Represents Relation Graph as defined in section 3.1 (Not generalised)_

## Constructor
`RelationGraph(l)`: The list of structures when the first element is the structure of N_0 and rest are block structures of N_i

## Fields
  * `vertices`: Set of vertices of the Relation Graph. The first element in the list is N_0 and rest are N_i.

  * `numOfVars`: Number of variables in the relation graph i.e. number of block structures.

## Functions
  * `isRelationGraph()`: Returns True iff the elements of the field vertices forms a relation graph.

  * `visualizeRelationGraph()`: Displays all the elements of the relation graph.

  * `findStrN0()`: Returns the compatible structure for N_0 given the set of vertices. Based on Algorithm 2.

  * `findStrN0Closure(U,s)`: Given a structure matrix U of size s*s it computes and returns the structure as described in algorithm 2.

  * `updateM(l,s,mat)`: Updation of the structure with the new found structure. Refer Algorithm 2.

  * `vertexCover(g,l,s)`: Returns the vertex cover of the graph G with incident edges on the set l and s as the number of nodes.

  * `calculateI(strN,size)`: Returns the set I for the structure strN of size ‘size’ as defined in Definition 3.2.1.

  * `findStrN0ForOne(blk)`: Implementation of Algorithm 1

  * `intersectionOfBlocks(X,Y,n)`: Returns the intersection of the two blocks X and Y of size n

# `Class Dynamics`
_Represents the dynamics of the system._

## Constructor
`Dynamics(A,n)`: The matrix A containing the dynamics of dimension n*n.

## Fields
  * `inpA`: contains the dynamics of the system i.e. A

  * `dimenstion`: A is of size dimension*dimension

## Functions
  * `displayA()`: Displays the matrix A.

  * `findConstStr()`: Returns the Structure (Definition 3.1) of N_0 (Definition 2.1) from the given dynamics.

  * `findConstMat()`: Returns N_0 (Definition 2.1) from the given dynamics.

  * `findVarBlkStr()`: Returns the block Structures of the variables in the system. As of now the paper doesn't explore which representation of the variable is the most optimal. Therefore, we just try to find a representation that doesn't take constant in any of the block structure. In some sense it is an exact one.
