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
import multiprocessing
import multiprocessing.pool
import time

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
        return (BlockStructure.intersection(self.getColStart(),self.getColBreadth(),mat.getRowStart(),mat.getRowLength())==False and BlockStructure.intersection(mat.getColStart(),mat.getColBreadth(),self.getRowStart(),self.getRowLength())==False)

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
        #print(self.vertices[1].mutualKerStrCheck(self.vertices[2]))
        #return True
        return self.checkStrN0()

    def checkStrN0(self):
        U=np.ones((self.vertices[1].size,self.vertices[1].size),dtype=int);

        for elem in self.vertices[1:]:
            U=RelationGraph.intersectionOfBlocks(U,RelationGraph.findStrN0ForOne(elem),self.vertices[1].size)

        for i in range(0,self.vertices[1].size):
            for j in range(0,self.vertices[1].size):
                if (U[i][j]==0 and self.vertices[0].matrix[i][j]!=0):
                    return False

        return RelationGraph.isStrClosed(U,self.vertices[1].size)

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
                print (self.inpA[i][j], " ", end=" ")
            print("\n")

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
        self.N0=np.zeros((self.dimension,self.dimension),dtype=float);
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
            ar=np.zeros((self.dimension,self.dimension),dtype=float)
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
        '''This is hardcoded for blocks of
        length and breadth 1 only'''
        M0=[]
        varMs=[]
        N0=self.matA.findConstMat()
        M0[:]=N0[:]
        varNs=self.matA.findVarBlkMats()
        varMBlkStrs=self.matA.findVarBlkStr()
        varMs[:]=varNs[:]
        for x in range(0,k-1):
            M0=np.matmul(N0,M0)
            for i in range(0,len(varMs)):
                c=varMBlkStrs[i].getRowStart()
                d=varMBlkStrs[i].getColStart()
                varMs[i][c][d]=(varNs[i][c][d]*M0[d][d])+(N0[c][c]*varMs[i][c][d])

        return [M0]+varMs

    def powerIntermidiate2(self, k):
        N0=M0=self.matA.findConstMat()
        varNs=varMs=self.matA.findVarBlkMats()
        for x in range(0,k-1):
            M0=np.matmul(N0,M0)
            for i in range(0,len(varMs)):
                varMs[i]=np.add(np.matmul(varNs[i],M0),np.matmul(N0,varMs[i]))

        return [M0]+varMs

    def power(self,k):
        A=np.zeros((self.matA.dimension,self.matA.dimension),dtype=float)
        for m in self.powerIntermidiate(k):
            A=np.add(A,m)
        return (Dynamics(A,self.matA.dimension,self.matA.vars))





class CalcReachableSet:

    def __init__(self,A,k,initset,flt,prtb):
        self.matA=A
        self.matAk=A.power(k)
        self.initialSet=initset
        self.faults=flt
        self.perturbations=prtb

    def formOpt(self):
        x=self.formX()
        pred=self.predicate()
        pert=self.perturbationsPredicate()

    def formX(self):
        V=np.zeros((self.matA.matA.dimension,self.matA.matA.dimension),dtype=float)
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

    class Stack:
     def __init__(self):
         self.items = []

     def isEmpty(self):
         return self.items == []

     def push(self, item):
         self.items.append(item)

     def pop(self):
         return self.items.pop()

     def peek(self):
         return self.items[len(self.items)-1]

     def size(self):
         return len(self.items)

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
        A=np.array([[4,5,0,0],[0,6,0,0],[0,8,1,2],[0,0,0,6]])
        A=np.array([[0,0,0,9],[2,6,0,0],[0,0,0,2],[0,0,0,6]])
        #A=np.array([[4,1,0,9],[6,3,0,0],[9,0,7,0],[7,0,9,8]])
        initset=[(1,2),(3,4),(5,6),(7,8)]
        for i in range(1,17):
            for elem in itertools.combinations(range(0,16),i):
                f=[]
                for e in elem:
                    f.append((e//4,e%4))
                inp=InputFormation(A,4,f,range(1,len(f)),initset,5)
                dynA=Dynamics(inp.formMatA(),4,f)
                rg=RelationGraph([dynA.findConstStr()]+dynA.findVarBlkStr())
                if (rg.isRelationGraph()==True):
                    print("Yo! ",i)
                    rg.visualizeRelationGraph()
                    #exit()
                    #print("mat[",e//4,"][",e%4,"]", end=", ")


    def experimentalGenerator4Help(par):
        A=par[0][0]
        n=par[0][1]
        r=par[1][0]
        c=par[1][1]
        print("Starting with (",r,",",c,")")
        stk=Examples.Stack()
        flt=[(r,c)]
        blocked=list(set(list(itertools.product([c],range(n)))+list(itertools.product(range(n),[r]))))
        stk.push((flt,blocked))
        count=1
        while (stk.isEmpty()==False):
            print("It: ",count)
            pe=stk.pop()
            flt=pe[0]
            blocked=list(set(pe[1]))
            inp=InputFormation(A,n,flt,range(1,len(flt)),range(1,len(flt)),5)
            dynA=Dynamics(inp.formMatA(),n,flt)
            rg=RelationGraph([dynA.findConstStr()]+dynA.findVarBlkStr())
            if (rg.isRelationGraph()):
                rg.visualizeRelationGraph()
            if (set(blocked)==set(itertools.product(range(n),range(n)))):
                print("PUTS")
                break
            possiblePlaces=[x for x in (list(itertools.product(range(n),range(n)))) if x not in blocked]

            for i in blocked:
                print(i)
            for index in possiblePlaces:
                blockedTmp=list(itertools.product([index[1]],range(n)))+list(itertools.product(range(n),[index[0]]))
                stk.push((flt+[index],blocked+blockedTmp))
            count+=1
        print("Finished (",r,",",c,")")
        print("\n\n")


    def experimentalGenerator4():
        #A=np.array([[4,0,0,9],[0,3,0,0],[9,0,7,0],[0,0,9,8]])
        #A=np.array([[4,1,0,9],[6,3,0,0],[9,0,7,0],[7,0,9,8]])
        A=np.array([[4,0,0,9],[0,3,0,0],[9,0,7,0],[0,0,9,8]])
        #initset=[(1,2),(3,4),(5,6),(7,8)]
        n=4
        '''inputrcs=list(itertools.product(range(n),range(n)))
        inputparams=list(itertools.product([(A,n)],inputrcs))
        pool = multiprocessing.Pool(8)
        pool.map(Examples.experimentalGenerator4Help,inputparams)'''


    def experimentalGenerator5Help(par):
        A=par[0][0]
        n=par[0][1]
        r=par[1][0]
        c=par[1][1]
        flt=[(r,c)]
        if A[r][c]!=0:
            inp=InputFormation(A,n,flt,range(1,len(flt)),range(1,len(flt)),5)
            dynA=Dynamics(inp.formMatA(),n,flt)
            rg=RelationGraph([dynA.findConstStr()]+dynA.findVarBlkStr())
            if (rg.isRelationGraph()):
                return (r,c)

    def experimentalGenerator5Help2(par):
        A=par[0][0]
        n=par[0][1]
        flts=par[0][2]
        nC=par[1]
        for elem in itertools.combinations(flts,nC):
            inp=InputFormation(A,n,elem,range(1,len(elem)),range(1,len(elem)),5)
            dynA=Dynamics(inp.formMatA(),n,elem)
            rg=RelationGraph([dynA.findConstStr()]+dynA.findVarBlkStr())
            if (rg.isRelationGraph()):
                return rg

    def experimentalGenerator5():
        #A=np.array([[0,0,0,9],[0,3,0,0],[0,0,7,0],[0,0,9,8]])
        #A=np.array([[0,0,0,9],[2,6,0,0],[0,0,0,2],[0,0,0,6]])
        '''A=np.array([[0,0,0,9,7,0],
                    [0,0,0,0,0,0],
                    [0,0,0,2,0,0],
                    [0,0,0,6,0,0],
                    [0,0,0,0,0,0],
                    [0,0,0,0,0,0]])'''
        #A=np.array([[5,0,8,9],[2,6,0,6],[9,0,8,2],[7,0,0,6]])
        '''A=np.array([[0,1,0,0,0,0,0,0,0,0,0,0,0,0,0],
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
                    ])'''
        '''A=np.array([[0,0,0,1,0,0,0,0,0,0,0,0],
                    [0,0,0,0,1,0,0,0,0,0,0,0],
                    [0,0,0,0,0,1,0,0,0,0,0,0],
                    [0,0,0,0,0,0,0,-9.8,0,0,0,0],
                    [0,0,0,0,0,0,9.8,0,0,0,0,0],
                    [0,0,0,0,0,0,0,0,0,0,0,0],
                    [0,0,0,0,0,0,0,0,0,1,0,0],
                    [0,0,0,0,0,0,0,0,0,0,1,0],
                    [0,0,0,0,0,0,0,0,0,0,0,1],
                    [0,0,0,0,0,0,0,0,0,0,0,0],
                    [0,0,0,0,0,0,0,0,0,0,0,0],
                    [0,0,0,0,0,0,0,0,0,0,0,0],
                    ])'''
        '''A=np.array([[0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0],
                    [0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0],
                    [0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0],
                    [0,0,0,0,0,0,0,-9.8,0,0,0,0,0,0,0,0],
                    [0,0,0,0,0,0,9.8,0,0,0,0,0,0,0,0,0],
                    [0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0],
                    [0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0],
                    [0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0],
                    [0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0],
                    [0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0],
                    [0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0],
                    [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1],
                    [0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0],
                    [0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0],
                    [0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0],
                    [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1],
                    ])'''
        '''A=np.array([[1,1,1,0],
                    [1,1,0,0],
                    [1,0,-1,0],
                    [1,0,0,-1],
                    ])'''
        '''A=np.array([[1,0],
                    [0,1],
                    ])'''
        '''A=np.array([[0,0,1,0],
                    [0,0,0,1],
                    [3.3,0,0,-9.3],
                    [0,0,-9,0],
                    ])'''
        '''A=np.array([[0,0,0],
                    [0,0,0],
                    [0,0,-1],
                    ])'''
        '''A=np.array([[-1,1,1,0,1],
                    [1,-1,0,0,0],
                    [1,0,-1,0,0],
                    [1,0,0,-1,0],
                    [0,0,0,0,1],
                    ])'''
        '''A=np.array([[0,0,1,0,0,0],
                    [0,0,0,1,0,0],
                    [1,0,0,1,1,0],
                    [0,0,-1,0,0,1],
                    [0,0,0,0,1,0],
                    [0,0,0,0,0,1],
                    ])'''
        '''A=np.array([[0,0,0,-1],
                    [0,0,0,0],
                    [0,0,-1,0],
                    [0,0,0,1],
                    ])'''
        '''A=np.array([[0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
                    [0,0,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
                    [1.702423734,3.929356551,-4.3607983776,-1.01374489,-1.6167727749,0.2653009364,-0.2375199245,-0.4793543458,0.06412815,-0.1079326841,-0.2463610381,0.0276872161,-0.0605561959,-0.1501445039,0.0151944922,-0.0374830081,-0.0986391305,0.009628751,-0.0242136837,-0.0665592904,0.0067836913,-0.015601062,-0.0442510048,0.0052325207,-0.0093924696,-0.0272127915,0.0043984935,-0.0044278796,-0.0129879863,0.0040303349],
                    [0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
                    [0,0,1,0,0,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
                    [0.688678844,2.3125837761,0.2653009364,1.4649038095,3.4500022052,-4.2966702275,-1.1216775741,-1.863133813,0.2929881525,-0.2980761204,-0.6294988497,0.0793226422,-0.1454156921,-0.3450001686,0.0373159671,-0.0847698796,-0.2167037943,0.0219781835,-0.0530840701,-0.1428901352,0.0148612718,-0.0336061533,-0.0937720819,0.0111821848,-0.0200289416,-0.057238991,0.0092628557,-0.0093924696,-0.0272127915,0.0084288284],
                    [0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
                    [0,0,0,0,0,1,0,0,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
                    [0.4511589195,1.8332294303,0.06412815,0.5807461599,2.066222738,0.2929881525,1.4043476136,3.2998577013,-4.2814757354,-1.1591605822,-1.9617729435,0.3026169036,-0.3222898041,-0.6960581401,0.0861063336,-0.1610167541,-0.3892511733,0.0425484878,-0.0941623492,-0.2439165858,0.026376677,-0.0575119497,-0.1558781215,0.0188916067,-0.0336061533,-0.0937720819,0.0152125197,-0.015601062,-0.0442510048,0.0136613491],
                    [0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
                    [0,0,0,0,0,0,0,0,1,0,0,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
                    [0.3432262354,1.5868683922,0.0276872161,0.3906027236,1.6830849264,0.0793226422,0.5432631518,1.9675836075,0.3026169036,1.3801339299,3.2332984109,-4.274692044,-1.1747616442,-2.0060239482,0.3078494243,-0.3316822737,-0.7232709316,0.090504827,-0.1654446337,-0.4022391596,0.0465788228,-0.0941623492,-0.2439165858,0.0304070119,-0.0530840701,-0.1428901352,0.0232901001,-0.0242136837,-0.0665592904,0.0204450405],
                    [0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
                    [0,0,0,0,0,0,0,0,0,0,0,1,0,0,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
                    [0.2826700395,1.4367238883,0.0151944922,0.3057432273,1.4882292617,0.0373159671,0.3663890398,1.616525636,0.0861063336,0.5276620899,1.9233326028,0.3078494243,1.3707414603,3.2060856194,-4.2702935506,-1.1791895238,-2.0190119345,0.3118797592,-0.3316822737,-0.7232709316,0.094535162,-0.1610167541,-0.3892511733,0.0509773162,-0.0847698796,-0.2167037943,0.0356395326,-0.0374830081,-0.0986391305,0.0300737915],
                    [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0],
                    [0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,-1,0,0,0,0,0,0,0,0,0,0,0,0],
                    [0.2451870315,1.3380847578,0.009628751,0.2584563558,1.3701645979,0.0219781835,0.2901421653,1.443978257,0.0425484878,0.3569965702,1.5893128445,0.090504827,0.5232342102,1.9103446165,0.3118797592,1.3707414603,3.2060856194,-4.2662632156,-1.1747616442,-2.0060239482,0.3162782527,-0.3222898041,-0.6960581401,0.0997676827,-0.1454156921,-0.3450001686,0.0577610076,-0.0605561959,-0.1501445039,0.0452682837],
                    [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0],
                    [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,-1,0,0,0,0,0,0,0,0,0],
                    [0.2209733477,1.2715254674,0.0067836913,0.2295859695,1.2938337531,0.0148612718,0.2490638862,1.3429518064,0.026376677,0.2857142857,1.4309902707,0.0465788228,0.3569965702,1.5893128445,0.094535162,0.5276620899,1.9233326028,0.3162782527,1.3801339299,3.2332984109,-4.2610306949,-1.1591605822,-1.9617729435,0.323061944,-0.2980761204,-0.6294988497,0.1093964337,-0.1079326841,-0.2463610381,0.0729554998],
                    [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0],
                    [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,-1,0,0,0,0,0,0],
                    [0.2053722857,1.2272744627,0.0052325207,0.2115808781,1.2443126759,0.0111821848,0.2251580898,1.2808457668,0.0188916067,0.2490638862,1.3429518064,0.0304070119,0.2901421653,1.443978257,0.0509773162,0.3663890398,1.616525636,0.0997676827,0.5432631518,1.9675836075,0.323061944,1.4043476136,3.2998577013,-4.2514019439,-1.1216775741,-1.863133813,0.3382564362,-0.2375199245,-0.4793543458,0.1370836498],
                    [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0],
                    [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,-1,0,0,0],
                    [0.1959798161,1.2000616712,0.0043984935,0.2009444061,1.2142864764,0.0092628557,0.2115808781,1.2443126759,0.0152125197,0.2295859695,1.2938337531,0.0232901001,0.2584563558,1.3701645979,0.0356395326,0.3057432273,1.4882292617,0.0577610076,0.3906027236,1.6830849264,0.1093964337,0.5807461599,2.066222738,0.3382564362,1.4649038095,3.4500022052,-4.2237147278,-1.01374489,-1.6167727749,0.4023845862],
                    [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0],
                    [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,-1],
                    [0.1915519365,1.1870736849,0.0040303349,0.1959798161,1.2000616712,0.0084288284,0.2053722857,1.2272744627,0.0136613491,0.2209733477,1.2715254674,0.0204450405,0.2451870315,1.3380847578,0.0300737915,0.2826700395,1.4367238883,0.0452682837,0.3432262354,1.5868683922,0.0729554998,0.4511589195,1.8332294303,0.1370836498,0.688678844,2.3125837761,0.4023845862,1.702423734,3.929356551,-3.9584137913]
                    ])'''
        '''A=np.array([[0,1,0,-1,0,0,-1,-1,0,0],
                    [0,0,1,0,0,0,0,0,0.3,0],
                    [0,0,-1,0,0,0,0,0,0.7,0],
                    [0,0,0,0,0,0.9,0,0,0,0],
                    [0,0,0,4,4,0,0.4,0.4,0,0],
                    [0,0,0,0,0,1,0,0,0,0],
                    [0,0,0,0,0,0,1,0,0,0],
                    [0,0,0,0,0,0,0,1,0,0],
                    [0,0,0,0,0,0,0,0,1,0],
                    [0,0,0,0,0,0,0,0,0,1],
                ])'''
        '''A=np.array([[0,0,0,0,0,1,0],
                    [0,0,0,0,0,0,1],
                    [1,0,0,0,0,0,0],
                    [0,1,0,0,0,0,0],
                    [0,0,0,0,0,0,0],
                    [0,0,0,0,0,1,0],
                    [0,0,0,0,0,0,1],
                    ])'''
        A=A=np.array([[0,1.0,0,-1.0,0,0,-0.00872664619,-0.00872664619,0,0],
                    [0,0,0.996,0,0,0,0,0,0.0003,0],
                    [0,0,-0.9999999525,0,0,0,0,0,0.07499999644,0],
                    [0,0,0,0,0,0.999999992,0,0,0,0],
                    [0,0,0,49.99999763,-49.99999763,0,0.4363322888,0.4363322888,0,0],
                    [0,0,0,0,0,1,0,0,0,0],
                    [0,0,0,0,0,0,1,0,0,0],
                    [0,0,0,0,0,0,0,1,0,0],
                    [0,0,0,0,0,0,0,0,1,0],
                    [0,0,0,0,0,0,0,0,0,1],
                    ])
        n=10
        inputrcs=list(itertools.product(range(n),range(n)))
        inputparams=list(itertools.product([(A,n)],inputrcs))
        pool = multiprocessing.Pool(8)
        l=pool.map(Examples.experimentalGenerator5Help,inputparams)
        l=list(filter(None, l))
        if not l:
            print("No fault possible")
            exit(0)
        print(l)
        print("----------")
        inputrcs=range(2,len(l)+1)
        inputparams=list(itertools.product([(A,n,l)],inputrcs))
        pool = multiprocessing.Pool(8)
        listRGs=pool.map(Examples.experimentalGenerator5Help2,inputparams)
        listRGs=list(filter(None, listRGs))
        for rg in listRGs:
            print("-----")
            rg.visualizeRelationGraph()
            print("-----")


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
        for i in range(1,15):
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

    def motorTransmission2():
        A=A=np.array([[0,0,0,0,0],
                    [7,0,0,0,0],
                    [0,0,1,0,0],
                    [0,0,0,1,0],
                    [3.2,3.2,0,0,1],
                    ])
        n=5
        f=[(4,0),(4,1)]
        v=[5,5] #yet to be decided
        T=20
        initset=[(0.9,1.1),(0.9,1.1),(0.9,1.1),(0.9,1.1),(0.9,1.1)] #Yet to figured out
        start_time = time.time()
        inpA=InputFormation(A,n,f,v,initset,T)
        dynA=Dynamics(inpA.formMatA(),n,f)
        rg=RelationGraph([dynA.findConstStr()]+dynA.findVarBlkStr())
        if (rg.isRelationGraph()==True):
            dynRepA=DynamicsRepOp(dynA)
            reach=CalcReachableSet(dynRepA,T,initset,f,inpA.perturbations)
            #----------------------
            end_time=time.time() - start_time
            print("---Motor Transmission Dynamics---")
            dynRepA.displayA()
            print("\n")
            print("---Matrix after time ",T,"---")
            reach.matAk.displayA()
            print("\n")
            print("Time taken: ",end_time,"\n")
        else:
            end_time=time.time() - start_time
            print("Out of scope of this Algorithm!")
            print("\nTime taken: ",end_time,"\n")

    def motorTransmission():
        A=A=np.array([[0,0,0,0,0],
                    [0,0,0,0,0],
                    [0,0,1,0,0],
                    [0,0,0,1,0],
                    [3.2,3.2,0,0,1],
                    ])
        n=5
        f=[(4,0),(4,1)]
        v=[5,5] #Yet to be decided
        T=20 #Yet to be figured out
        initset=[0,0,-0.0165,0.003,0] #Yet to figured out, says uncertain
        Examples.exampleTemplate(A,n,f,v,T,initset,"Motor Transmission")

    def holes():
        A=A=np.array([[0,1.0,0,-1.0,0,0,-0.00872664619,-0.00872664619,0,0],
                    [0,0,0.996,0,0,0,0,0,0.0003,0],
                    [0,0,-0.9999999525,0,0,0,0,0,0.07499999644,0],
                    [0,0,0,0,0,0.999999992,0,0,0,0],
                    [0,0,0,49.99999763,-49.99999763,0,0.4363322888,0.4363322888,0,0],
                    [0,0,0,0,0,1,0,0,0,0],
                    [0,0,0,0,0,0,1,0,0,0],
                    [0,0,0,0,0,0,0,1,0,0],
                    [0,0,0,0,0,0,0,0,1,0],
                    [0,0,0,0,0,0,0,0,0,1],
                    ])
        n=10
        f=[(0,6),(0,7),(4,6),(4,7)]
        v=[5,5,6,7] #Yet to be decided
        T=5 #Yet to be figured out
        initset=[(0.9,1.1),(0.9,1.1),(0.9,1.1),(0.9,1.1),(0.9,1.1),(0.9,1.1),(0.9,1.1),(0.9,1.1),(0.9,1.1),(0.9,1.1)] #Yet to figured out
        Examples.exampleTemplate(A,n,f,v,T,initset,"Holes Controller - XC")

    def transformer():
        A=A=np.array([[0,0,0,-9.63387],
                    [0,0,0,0],
                    [0,0,-36.354,0],
                    [0,0,0,1],
                    ]) #Yet to figured out
        n=4
        f=[(0,3)]
        v=[5] #Yet to be decided
        T=20 #Yet to be figured out
        initset=[(0.9,1.1),(0.9,1.1),(0.9,1.1),(0.9,1.1),(0.9,1.1)] #Yet to figured out
        Examples.exampleTemplate(A,n,f,v,T,initset,"Transformer Isolated DC-DC Converter")


    def exampleTemplate(A,n,f,v,T,initset,strName):
        start_time = time.time()
        inpA=InputFormation(A,n,f,v,initset,T)
        dynA=Dynamics(inpA.formMatA(),n,f)
        rg=RelationGraph([dynA.findConstStr()]+dynA.findVarBlkStr())
        if (rg.isRelationGraph()==True):
            dynRepA=DynamicsRepOp(dynA)
            reach=CalcReachableSet(dynRepA,T,initset,f,inpA.perturbations)
            #----------------------
            end_time=time.time() - start_time
            print("----",strName," Dynamics ----")
            dynRepA.displayA()
            print("\n")
            print("---- Matrix after time ",T," ----")
            reach.matAk.displayA()
            print("\n")
            print("Time taken: ",end_time,"\n")
        else:
            end_time=time.time() - start_time
            print("Out of scope of this Algorithm!")
            print("\nTime taken: ",end_time,"\n")




#Examples.experimentalGenerator5()
#Examples.motorTransmission()
Examples.holes()
#Examples.transformer()
'''A=np.array([[4,5,0,9],[0,6,0,0],[0,8,1,2],[0,0,0,6]])
f=[(2,3)]
initset=[(1,2),(3,4),(5,6),(7,8)]
inp=InputFormation(A,4,f,range(1,len(f)),initset,5)
dynA=Dynamics(inp.formMatA(),4,f)
rg=RelationGraph([dynA.findConstStr()]+dynA.findVarBlkStr())
if (rg.isRelationGraph()):
    rg.visualizeRelationGraph()
rg.findStrN0()'''

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
