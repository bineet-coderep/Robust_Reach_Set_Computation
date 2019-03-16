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
        flag=True
        for mat in self.vertices[1:]:
            for mat2 in self.vertices[1:]:
                #if (mat!=mat2):
                if mat.mutualKerStrCheck(mat2) == False:
                    return False
        return self.checkStrN0()

    def checkStrN0(self):
        U=np.ones((self.vertices[1].size,self.vertices[1].size),dtype=int);

        for elem in self.vertices[1:]:
            U=RelationGraph.intersectionOfBlocks(U,RelationGraph.findStrN0ForOne(elem),self.vertices[1].size)

        if np.array_equal(U,self.vertices[0].matrix):
            return RelationGraph.isStrClosed(U,self.vertices[1].size)
        return False

    @staticmethod
    def isStrClosed(U,n):
        U2=np.matmul(U,U)
        for i in range(0,n):
            for j in range(0,n):
                if U[i][j]==0 and U2[i][j]!=0:
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
                print (U[i][j], "  ", end="")
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

    def __init__(self, A, n, v):
        self.inpA=A; #contains the dynamics of the system i.e. A
        self.dimension=n; # A is of size dimension*dimension
        self.vars=v; #vars are the indices where there are perturbations

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
                if ((i,j) in self.vars or self.inpA[i][j]==0):
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
        '''
        self.N0=np.zeros((self.dimension,self.dimension),dtype=int);
        for i in range(0,self.dimension):
            for j in range(0,self.dimension):
                if ((i,j) not in self.vars):
                    self.N0[i][j]=self.inpA[i][j]
        return self.N0

    def findVarBlkStr(self):
        '''Finds out the block Structures of
        the variables in the system.
        Return Type: List of BlockStructure.
        '''
        l=[]
        for v in self.vars:
            l.append(BlockStructure((v[0],1,v[1],1),self.dimension))
        return l

    def findVarBlkMats(self):
        '''Finds the block structure matrices'''
        l=[]
        for v in self.findVarBlkStr():
            ar=np.zeros((self.dimension,self.dimension),dtype=int)
            ar[v.getRowStart()][v.getColStart()]=self.inpA[v.getRowStart()][v.getColStart()]
            l.append(ar)
        return l




class DynamicsRepOp:
    '''Represents the matrix as defined in Defn 2.1
    and performs operation as defined in Defn 2.1
    '''

    def __init__(self,A):
        self.matA=A;

    def displayA(self):
        N0=self.matA.findConstMat()
        varNs=self.matA.findVarBlkMats()
        for i in range(0,self.matA.dimension):
            for j in range(0,self.matA.dimension):
                print (N0[i][j],"  ",end="")
            print("\n")

        print("\n")
        for var in varNs:
            for i in range(0,self.matA.dimension):
                for j in range(0,self.matA.dimension):
                    print (var[i][j],"  ",end="")
                print("\n")
            print("\n")

    def powerIntermidiate(self, k):
        N0=M0=self.matA.findConstMat()
        varNs=varMs=self.matA.findVarBlkMats()
        for x in range(0,k-1):
            M0=np.matmul(N0,M0)
            for i in range(0,len(varMs)):
                varMs[i]=np.add(np.matmul(varNs[i],M0),np.matmul(N0,varMs[i]))

        return [M0]+varMs

    def power(self,k):
        A=np.zeros((self.matA.dimension,self.matA.dimension),dtype=int)
        for m in self.powerIntermidiate(k):
            A=np.add(A,m)
        return (Dynamics(A,self.matA.dimension,self.matA.vars))





class CalcReachableSet:

    def __init__(self,A,k,initset,prtb):
        self.matA=A
        self.matAk=A.power(k)
        self.initialSet=initset
        self.perturbations=prtb

    def formOpt(self):
        x=self.formX()
        pred=self.predicate()
        pert=self.perturbationsPredicate()

    def formX(self):
        V=np.zeros((self.matA.matA.dimension,self.matA.matA.dimension),dtype=int)
        for i in range(0,self.matA.matA.dimension):
            V[i][i]=1;
        X=np.matmul(self.matAk.inpA,V)
        return X

    def predicate(self):
        return self.initialSet

    def perturbationsPredicate(self):
        return self.perturbations;



class InputFormation:

    def __init__(self,A, n, v, vrnc, initset, k):
        self.dynA=A
        self.dimension=n
        self.vars=v
        self.variance=vrnc
        self.initialSet=initset
        self.time=k

    def formMatA(self):
        A=np.copy(self.dynA)
        for i in self.vars:
            A[i[0]][i[1]]=1

        return A


    def perturbations(self):
        pr=[]
        for i in range(0,len(self.vars)):
            t=(self.dynA[self.vars[i][0]][self.vars[i][1]]*self.variance[i])/100
            t2=self.dynA[self.vars[i][0]][self.vars[i][1]]+t
            t3=self.dynA[self.vars[i][0]][self.vars[i][1]]-t
            pr.append((t3,t2))

        return pr



class Examples:

    @staticmethod
    def experimental():
        #inpA=np.array([[2,0,9999,9999],[0,3,9999,9999],[0,0,6,7],[0,0,8,9]])
        A=np.array([[4,0,0,9],[0,3,0,0],[9,0,7,0],[0,0,9,8]])
        #f=[(0,3),(2,0),(3,1),(3,2)]
        f=[(0,3)]
        v=[5,5,5,5]
        initset=[(1,2),(3,4),(5,6),(7,8)]
        inp=InputFormation(A,4,f,v,initset,5)
        dynA=Dynamics(inp.formMatA(),4,f)
        rg=RelationGraph([dynA.findConstStr()]+dynA.findVarBlkStr())
        rg.visualizeRelationGraph()
        print(rg.isRelationGraph())
        if (rg.isRelationGraph()):
            ARep=DynamicsRepOp(dynA)
            rs=CalcReachableSet(ARep,2,initset,inp.perturbations())
            rs.formOpt()

    def experimentalGenerator():
        for i in range(1,16):
            for elem in itertools.combinations(range(1,226),i):
                print(len(elem))
                for e in elem:
                    l=i
                    #print("mat[",e//4,"][",e%4,"]", end=", ")

    def experimentalGenerator2():
        A=np.array([[4,0,0,9],[0,3,0,0],[9,0,7,0],[0,0,9,8]])
        initset=[(1,2),(3,4),(5,6),(7,8)]
        for i in range(4,17):
            for elem in itertools.combinations(range(0,16),i):
                f=[]
                for e in elem:
                    f.append((e//4,e%4))
                inp=InputFormation(A,4,f,range(1,len(f)),initset,5)
                dynA=Dynamics(inp.formMatA(),4,f)
                rg=RelationGraph([dynA.findConstStr()]+dynA.findVarBlkStr())
                if (rg.isRelationGraph()):
                    print("Yo! ",i)
                    rg.visualizeRelationGraph()
                    exit()
                    #print("mat[",e//4,"][",e%4,"]", end=", ")

    def experimentalGenerator3():
        A=np.array([[0,1,0,0,0,0,0,0,0,0,0,0,0,0,0],
                    [0,0,-1,0,0,0,0,0,0,0,0,0,0,0,0],
                    [1.7152555329,3.9705119979,-4.3600526739,-0.9999330812,-1.5731541104,0.2669165553,-0.2215507198,-0.4303855023,0.0669078193,-0.0881500219,-0.1881468451,0.0322187056,-0.0343095071,-0.0767587194,0.0226660281],
                    [0,0,0,0,1,0,0,0,0,0,0,0,0,0,0],
                    [0,0,1,0,0,-1,0,0,0,0,0,0,0,0,0],
                    [0.7153224517,2.3973578876,0.2669165553,1.4937048131,3.5401264957,-4.2931448546,-1.0880831031,-1.7613009555,0.2991352608,-0.2558602268,-0.5071442217,0.0895738474,-0.0881500219,-0.1881468451,0.0548847337],
                    [0,0,0,0,0,0,0,1,0,0,0,0,0,0,0],
                    [0,0,0,0,0,1,0,0,-1,0,0,0,0,0,0],
                    [0.493771732,1.9669723853,0.0669078193,0.6271724298,2.2092110425,0.2991352608,1.4593953061,3.4633677762,-4.2704788265,-1.0880831031,-1.7613009555,0.3218012889,-0.2215507198,-0.4303855023,0.121792553],
                    [0,0,0,0,0,0,0,0,0,0,1,0,0,0,0],
                    [0,0,0,0,0,0,0,0,1,0,0,-1,0,0,0],
                    [0.40562171,1.7788255402,0.0322187056,0.4594622249,1.8902136659,0.0895738474,0.6271724298,2.2092110425,0.3218012889,1.4937048131,3.5401264957,-4.2382601209,-0.9999330812,-1.5731541104,0.3887091083],
                    [0,0,0,0,0,0,0,0,0,0,0,0,0,1,0],
                    [0,0,0,0,0,0,0,0,0,0,0,1,0,0,-1],
                    [0.371312203,1.7020668208,0.0226660281,0.40562171,1.7788255402,0.0548847337,0.493771732,1.9669723853,0.121792553,0.7153224517,2.3973578876,0.3887091083,1.7152555329,3.9705119979,-3.9713435656]
                    ])
        initset=[(0.9,1.1),(0.9,1.1),(0.9,1.1),(0.9,1.1),(0.9,1.1),(0.9,1.1),(0.9,1.1),(0.9,1.1),(0.9,1.1),(0.9,1.1),(0.9,1.1),(0.9,1.1),(0.9,1.1),(0.9,1.1),(0.9,1.1)]
        for i in range(6,10):
            print("i: ",i)
            for elem in itertools.combinations(range(0,225),i):
                f=[]
                for e in elem:
                    f.append((e//15,e%15))
                inp=InputFormation(A,15,f,range(1,len(f)),initset,5)
                dynA=Dynamics(inp.formMatA(),15,f)
                rg=RelationGraph([dynA.findConstStr()]+dynA.findVarBlkStr())
                if (rg.isRelationGraph()):
                    print("Yo! ",i)
                    rg.visualizeRelationGraph()
                    exit()
                    #print("mat[",e//4,"][",e%4,"]", end=", ")


    @staticmethod
    def fiveVehiclePlatoon():
        A=np.array([[0,1,0,0,0,0,0,0,0,0,0,0,0,0,0],
                    [0,0,-1,0,0,0,0,0,0,0,0,0,0,0,0],
                    [1.7152555329,3.9705119979,-4.3600526739,-0.9999330812,-1.5731541104,0.2669165553,-0.2215507198,-0.4303855023,0.0669078193,-0.0881500219,-0.1881468451,0.0322187056,-0.0343095071,-0.0767587194,0.0226660281],
                    [0,0,0,0,1,0,0,0,0,0,0,0,0,0,0],
                    [0,0,1,0,0,-1,0,0,0,0,0,0,0,0,0],
                    [0.7153224517,2.3973578876,0.2669165553,1.4937048131,3.5401264957,-4.2931448546,-1.0880831031,-1.7613009555,0.2991352608,-0.2558602268,-0.5071442217,0.0895738474,-0.0881500219,-0.1881468451,0.0548847337],
                    [0,0,0,0,0,0,0,1,0,0,0,0,0,0,0],
                    [0,0,0,0,0,1,0,0,-1,0,0,0,0,0,0],
                    [0.493771732,1.9669723853,0.0669078193,0.6271724298,2.2092110425,0.2991352608,1.4593953061,3.4633677762,-4.2704788265,-1.0880831031,-1.7613009555,0.3218012889,-0.2215507198,-0.4303855023,0.121792553],
                    [0,0,0,0,0,0,0,0,0,0,1,0,0,0,0],
                    [0,0,0,0,0,0,0,0,1,0,0,-1,0,0,0],
                    [0.40562171,1.7788255402,0.0322187056,0.4594622249,1.8902136659,0.0895738474,0.6271724298,2.2092110425,0.3218012889,1.4937048131,3.5401264957,-4.2382601209,-0.9999330812,-1.5731541104,0.3887091083],
                    [0,0,0,0,0,0,0,0,0,0,0,0,0,1,0],
                    [0,0,0,0,0,0,0,0,0,0,0,1,0,0,-1],
                    [0.371312203,1.7020668208,0.0226660281,0.40562171,1.7788255402,0.0548847337,0.493771732,1.9669723853,0.121792553,0.7153224517,2.3973578876,0.3887091083,1.7152555329,3.9705119979,-3.9713435656]
                    ])
        f=[(3,2)]
        v=[5]
        initset=[(0.9,1.1),(0.9,1.1),(0.9,1.1),(0.9,1.1),(0.9,1.1),(0.9,1.1),(0.9,1.1),(0.9,1.1),(0.9,1.1),(0.9,1.1),(0.9,1.1),(0.9,1.1),(0.9,1.1),(0.9,1.1),(0.9,1.1)]
        inpA=InputFormation(A,15,f,v,initset,5)
        dynA=Dynamics(inpA.formMatA(),15,f)
        rg=RelationGraph([dynA.findConstStr()]+dynA.findVarBlkStr())
        rg.visualizeRelationGraph()
        print("---")
        rg.findStrN0()




Examples.experimentalGenerator3();
#inpA=np.array([[2,0,9999,9999],[0,3,9999,9999],[0,0,6,7],[0,0,8,9]])
'''inpA=np.array([[4,0,0,9],[0,3,0,0],[9,0,7,0],[0,0,9,8]])
A=Dynamics(inpA,4,[(0,3),(2,0),(3,1),(3,2)],5)
ARep=DynamicsRepOp(A)
ARep.displayA()
Ak=ARep.power(3)
Ak.displayA()
print();
#A.displayA()
#print()
rg=RelationGraph([A.findConstStr()]+A.findVarBlkStr())
#rg.visualizeRelationGraph()
#print(rg.isRelationGraph())
#print()
ll=rg.findStrN0()'''
'''for i in range(0,4):
    for j in range(0,4):
        print (ll[i][j], "      ", end="");
    print();'''

#print()


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
