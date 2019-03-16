'''
Author: Bineet Ghosh, under supervision of Prof. Sridhar
Email: ghosh.bineet22@gmail.com
Date: April 19,2018

Based On Reachable set computation of
Linear Systems with Faults using
Block Structure Matrices

API Documentation: API_Doc.md
'''

import numpy as np
import itertools

class Structure:
    '''Represents structure of a matrix
    Defn 3.1
    '''

    #As of now only the non-zero tuples are stored, with no provision of multiple sets
    def __init__(self,l,n):
        self.omegas=l #List uple of indices: (startrow,length,startcol,breadth)
        self.size=n
        '''creates an actual matrix Representing this structure
        1 represents the set R, and 0 otherwise '''
        self.matrix=np.zeros((self.size,self.size),dtype=int)
        for omega in self.omegas:
            for i in range(omega[0],omega[0]+omega[1]):
                for j in range(omega[2],omega[2]+omega[3]):
                    self.matrix[i][j]=1 #1 Represents the set R

    def visualizeStr(self):
        for i in range(0,self.size):
            for j in range(0,self.size):
                print (self.matrix[i][j], " ", end="")
            print ("\n")


class BlockStructure(Structure):
    '''Represents block structured Matrices
    Defn 3.3
    '''

    def __init__(self,l,n):
        Structure.__init__(self,[l],n);

    def getRowStart(self):
        return self.omegas[0][0]

    def getRowLength(self):
        return self.omegas[0][1]

    def getColStart(self):
        return self.omegas[0][2]

    def getColBreadth(self):
        return self.omegas[0][3]

    @staticmethod
    def intersection(x,lx,y,ly):
        #[[x,lx]] \cap [[y,ly]]
        if (x<=y):
            return (x+lx-1 >= y)
        else:
            return (y+ly-1 >= x)

    def mutualKerStrCheck(self,mat):
        #Defn 3.4
        return not (BlockStructure.intersection(self.getColStart(),self.getColBreadth(),mat.getRowStart(),mat.getRowLength()) and BlockStructure.intersection(mat.getColStart(),mat.getColBreadth(),self.getRowStart(),self.getRowLength()))

    #def generateAll()

class RelationGraph:
    '''Represents the relation graph
    as defined in section 3.1 (Not genralized)
    '''
    def __init__(self, l):
        self.vertices=l #l[0]=N_0 and rest l[i]=N_i
        self.numOfVars=len(l)-1

    def isRelationGraph(self):
        #Here self.vertices are <N> and <M> as well.
        #N_i \times M_i = 0: Not required as of now with the current setup
        #Checks if (N_i \times M_j) + (M_i \times N_j) = 0
        for mat in self.vertices[1:]:
            for mat2 in self.vertices[1:]:
                #if (mat!=mat2):
                if mat.mutualKerStrCheck(mat2) == False:
                    return False
        return True

    def visualizeRelationGraph(self):
        self.vertices[0].visualizeStr() #Structure of N_0
        print()
        for elem in self.vertices[1:]:
            elem.visualizeStr() #BlockStructure of N_i
            print()

    def findStrN0(self):
        U=np.ones((self.vertices[1].size,self.vertices[1].size),dtype=int);

        for elem in self.vertices[1:]:
            U=RelationGraph.intersectionOfBlocks(U,RelationGraph.findStrN0ForOne(elem),self.vertices[1].size)
            #U=Dynamics.findStrN0ForOne(elem);

        RelationGraph.findStrN0Closure(U,self.vertices[1].size)

        for i in range(0,self.vertices[1].size):
            for j in range(0,self.vertices[1].size):
                print (U[i][j], "      ", end="")
            print()

    def findStrN0Closure(U,s):
        #Algorithm 2
        varI=RelationGraph.calculateI(U,s)
        M=np.copy(U)
        while (varI):
            graphG=np.zeros((s*s,s*s),dtype=int) #1-D representation of a 2-D array is used
            ver=[] #Set of participating vertices;
            for c in varI:
                for i in range(0,s):
                    if (U[c[0]][i]!=0 and U[i][c[1]]!=0):
                        graphG[c[0]*s+i][i*s+c[1]]=1 #1-D representation of a 2-D array is used
                        ver.append(c[0]*s+i)
                        ver.append(i*s+c[1])
            vc=RelationGraph.vertexCover(graphG,ver,s*s)
            M=RelationGraph.updateM(vc,s,M)
            varI=RelationGraph.calculateI(M,s)

    def updateM(l,s,mat):
        '''Updates the structure matrix as defined
        in Algorithm 2'''
        for i in l:
            mat[i//s][i%s]=0
        return mat

    def vertexCover(g,l,s):
        for i in range(1,s+1):
            for elem in itertools.combinations(l,i):
                #is 'elem' vertex cover of g?
                gtemp=np.copy(g);
                for v in elem:
                    for j in range(0,s):
                        gtemp[j][v]=0
                        gtemp[v][j]=0
                flag=True;
                for j in range(0,s):
                    for k in range(0,s):
                        if gtemp[j][k]==1:
                            flag=False
                            break
                if flag:
                    return elem;

                print;

    def calculateI(strN,size):
        '''Calculates the set I according to
        definition 3.2.1'''
        t=np.matmul(strN,strN)
        l=[] #Set I as defined in 3.2.1
        for i in range(0,size):
            for j in range(0,size):
                if (strN[i][j]==0 and t[i][j]!=0):
                    l.append((i,j))
        return l


    @staticmethod
    def findStrN0ForOne(blk):
        #Algorithm 1
        M=np.ones((blk.size,blk.size),dtype=int)
        O=np.ones((blk.size,blk.size),dtype=int)

        for i in range(blk.getRowStart(),blk.getRowStart()+blk.getRowLength()):
            for j in range(0,blk.size):
                O[j][i]=0

        for i in range(blk.getRowStart(),blk.getRowStart()+blk.getRowLength()):
            for j in range(0,blk.size):
                O[i][j]=1

        for i in range(blk.getColStart(),blk.getColStart()+blk.getColBreadth()):
            for j in range(0,blk.size):
                M[i][j]=0

        for i in range(blk.getColStart(),blk.getColStart()+blk.getColBreadth()):
            for j in range(0,blk.size):
                M[j][i]=1

        return (RelationGraph.intersectionOfBlocks(O,M,blk.size))

    @staticmethod
    def intersectionOfBlocks(X,Y,n):
        #Calculates intersection of two blocks
        Z=np.empty((n,n),dtype=int);
        for i in range(0,n):
            for j in range(0,n):
                Z[i][j]=X[i][j]*Y[i][j]
        return Z;


class Dynamics:
    '''Represnts the dynamics of the System'''


    def __init__(self, A, n):
        self.inpA=A; #contains the dynamics of the system i.e. A
        self.dimension=n; # A is of size dimension*dimension

    def displayA(self):
        for i in range(0,self.dimension):
            for j in range(0,self.dimension):
                print (self.inpA[i][j], "      ", end="")
            print()

    def findConstStr(self):
        '''Finds out the Structure (Definition 3.1)
        of N_0 (Definition 2.1) from the given dynamics.
        Return Type: Structure
        '''
        l=[]
        for i in range(0,self.dimension):
            startrow=i
            startcol=0
            breadth=0
            for j in range(0,self.dimension):
                if (inpA[i][j]==9999 or inpA[i][j]==0):
                    if (breadth>0):
                        l.append((startrow,1,startcol,breadth))
                    breadth=0
                    startcol=j+1
                else:
                    breadth+=1
            if (breadth>0):
                l.append((startrow,1,startcol,breadth))

        return (Structure(l,self.dimension))

    def findConstMat(self):
        '''Finds out the N_0 (Definition 2.1)
        from the given dynamics.
        Return Type: Structure
        '''
        self.N0=np.zeros((self.dimension,self.dimension));
        for i in range(0,self.dimension):
            for j in range(0,self.dimension):
                if (inpA[i][j]!=9999):
                    self.N0[i][j]=inpA[i][j]

    def findVarBlkStr(self):
        '''Finds out the block Structures of
        the variables in the system.
        Return Type: List of BlockStructure.
        Implementation: As of now the paper doesn't explore
        which representation of the variable is the most
        optimal. Therefore, we just try to find a representation
        that doesn't take constant in any of the block structure.
        In some sense it is an exact one.
        Note to Self/Author: Exploration reqd, as this could be the optimal
        (At least one of them)
        '''
        l=[]
        for i in range(0,self.dimension):
            startrow=i
            startcol=0
            breadth=0
            for j in range(0,self.dimension):
                if (inpA[i][j]!=9999):
                    if (breadth>0):
                        l.append(BlockStructure((startrow,1,startcol,breadth),self.dimension))
                    breadth=0
                    startcol=j+1
                else:
                    breadth+=1
            if (breadth>0):
                l.append(BlockStructure((startrow,1,startcol,breadth),self.dimension))
        return l







#inpA=np.array([[2,0,9999,9999],[0,3,9999,9999],[0,0,6,7],[0,0,8,9]])
inpA=np.array([[4,3,4,9999],[0,3,6,4],[9999,0,6,8],[6,9999,9999,9]])
A=Dynamics(inpA,4)
A.displayA()
print()
rg=RelationGraph([A.findConstStr()]+A.findVarBlkStr())
rg.visualizeRelationGraph()
print(rg.isRelationGraph())
print()
ll=rg.findStrN0()
'''for i in range(0,4):
    for j in range(0,4):
        print (ll[i][j], "      ", end="");
    print();'''

print()


'''a=Structure([(0,1,3,1),(2,2,3,1)],5);
a.visualizeMat();
print();
b=BlockStructure((0,3,1,3),5);
b.visualizeMat();
print("");
c=BlockStructure((0,3,1,3),5);
c.visualizeMat();

g=RelationGraph([a,b,c],2);
print(g.isRelationGraph());
'''
