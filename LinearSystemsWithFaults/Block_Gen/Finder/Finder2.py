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
import math
from gurobipy import *

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

    def getFaultyCells(self):
        l=[]
        for i in range(self.getRowStart(),self.getRowStart()+self.getRowLength()):
            for j in range(self.getColStart(),self.getColStart()+self.getColBreadth()):
                l.append((i,j))
        return l


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
        U=np.ones((self.vertices[1].size,self.vertices[1].size),dtype=int)

        for elem in self.vertices[1:]:
            U=RelationGraph.intersectionOfBlocks(U,RelationGraph.findStrN0ForOne(elem),self.vertices[1].size)

        for i in range(0,self.vertices[1].size):
            for j in range(0,self.vertices[1].size):
                if (U[i][j]==0 and self.vertices[0].matrix[i][j]!=0):
                    return False

        #print(U)
        return RelationGraph.isStrClosed(U,self.vertices[1].size)
        #return RelationGraph.isStrClosed(self.vertices[0].matrix,self.vertices[0].size)

    @staticmethod
    def isStrClosed(U,n):
        U2=np.matmul(U,U)
        for i in range(0,n):
            for j in range(0,n):
                if U[i][j]==0 and U2[i][j]!=0:
                    return False

        #print(U)
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
    '''Represents the dynamics of the System'''

    def __init__(self, A, n, v):
        self.inpA=A; #contains the dynamics of the system i.e. A
        self.dimension=n; # A is of size dimension*dimension
        self.vars=v; #set of blocks indicating faults

    def displayA(self):
        for i in range(0,self.dimension):
            for j in range(0,self.dimension):
                print (self.inpA[i][j], " ", end=" ")
            print("\n")

    def findAllFaultyCells(self):
        b=self.findVarBlkStr()
        fl=[]
        for e in b:
            fl=fl+e.getFaultyCells()

        return fl

    def findConstStr(self):
        '''Finds out the Structure (Definition 3.1)
        of N_0 (Definition 2.1) from the given dynamics.
        Return Type: Structure
        '''
        fltList=self.findAllFaultyCells()
        l=[]
        for i in range(0,self.dimension):
            startrow=i
            startcol=0
            breadth=0
            for j in range(0,self.dimension):
                if (((i,j) in fltList) or self.inpA[i][j]==0):
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
        fltList=self.findAllFaultyCells()
        self.N0=np.zeros((self.dimension,self.dimension),dtype=np.float128);
        for i in range(0,self.dimension):
            for j in range(0,self.dimension):
                if ((i,j) not in fltList):
                    self.N0[i][j]=self.inpA[i][j]
        return self.N0

    def findVarBlkStr(self):
        '''Finds out the block Structures of
        the variables in the system.
        Return Type: List of BlockStructure.
        '''
        l=[]
        for v in self.vars:
            l.append(BlockStructure(v,self.dimension))
        return l


    def findVarBlkMats(self):
        '''Finds the block structure matrices'''
        l=[]
        for v in self.findVarBlkStr():
            ar=np.zeros((self.dimension,self.dimension),dtype=np.float128)
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
                varMs[i]=np.matmul(varMs[i],varNs[i])


        return [M0]+varMs

    def power(self,k):
        A=np.zeros((self.matA.dimension,self.matA.dimension),dtype=np.float128)
        for m in self.powerIntermidiate(k):
            A=np.add(A,m)
        return (Dynamics(A,self.matA.dimension,self.matA.vars))

class InputFormation:

    def __init__(self,A, n, v, vrnc):
        self.dynA=A
        self.dimension=n
        self.vars=v
        self.variance=vrnc

    def formMatA(self):
        f=self.getAllFaultyCells()
        A=np.copy(self.dynA)
        for i in f:
            A[i[0]][i[1]]=1

        return A

    def getAllFaultyCells(self):
        l=[]
        for e in self.vars:
            for i in range(e[0],e[0]+e[1]):
                for j in range(e[2],e[2]+e[3]):
                    l.append((i,j))
        return l

    def perturbations(self):
        pr=[]
        for i in range(len(self.vars)):
            mn=minimumInCells(self,i)
            mx=maximumInCells(self,i)
            t1=mn-((mn*self.variance[i])/100)
            t2=mx+((mx*self.variance[i])/100)
            pr.append(min(t1,t2),max(t1,t2))

        return pr

    def minimumInCells(self,i):
        mn=9999999
        for i in rang(e[0],e[0]+e[1]):
            for j in range(e[2],e[2]+e[3]):
                if self.dynA[i][j]<mn:
                    mn=self.dynA[i][j]
        return mn

    def maximumInCells(self,i):
        mx=-9999999
        for i in rang(e[0],e[0]+e[1]):
            for j in range(e[2],e[2]+e[3]):
                if self.dynA[i][j]>mx:
                    mx=self.dynA[i][j]
        return mx

class CalcReachableSet:

    def __init__(self,A,k,initset,flt,prtb):
        self.matA=A
        self.matAk=A.power(k)
        self.initialSet=initset
        self.faults=flt
        self.perturbations=prtb

    def getAllFaultyCells(self):
        l=[]
        for e in self.faults:
            for i in rang(e[0],e[0]+e[1]):
                for j in range(e[2],e[2]+e[3]):
                    l.append((i,j))
        return l

    def indexOfBlk(self, i,j):
        for ind in range(len(self.faults)):
            if i>=e[0] and i<e[0]+e[1]:
                if j>=e[2] and j<e[2][3]:
                    return ind

    def isFeasible(self,A,B,ln):
        flag=True
        model = Model("qp")
        #model.Params.BarHomogeneous=1
        #model.Params.NumericFocus=3

        #create variables
        faultVars=[]
        for i in self.faults:
            name="fault "+str(i)
            faultVars.append(model.addVar(-GRB.INFINITY,GRB.INFINITY,name=name,vtype='C'))
        dynVars=[]
        for i in range(self.matA.matA.dimension):
            name="alpha"+str(i)
            dynVars.append(model.addVar(-GRB.INFINITY,GRB.INFINITY,name=name,vtype='C'))

        #Add initial Set and Perturbations constraints
        for i in range(self.matA.matA.dimension):
            model.optimize()
            name="Init-C"+str(i)
            model.addConstr(dynVars[i]>=self.initialSet[i][0],name+".1")
            model.addConstr(dynVars[i]<=self.initialSet[i][1],name+".2")
        no=0
        for i in faultVars:
            name="Fault-C"+str(no)
            print(self.perturbations[no][0],",",self.perturbations[no][1])
            model.addConstr(i>=self.perturbations[no][0],name+".1")
            model.addConstr(i<=self.perturbations[no][1],name+".2")
            no+=1

        #Generate Final Equations
        rVars=[]
        for i in range(self.matA.matA.dimension):
            obj=0
            for j in range(self.matA.matA.dimension):
                if (i,j) in self.getAllFaultyCells():
                    obj=obj+(self.matAk.inpA[i][j]*faultVars[self.indexOfBlk(i,j)]*dynVars[j])
                else:
                    obj=obj+(self.matAk.inpA[i][j]*dynVars[j])
            rVars.append(obj)


        '''for i in range(ln):
            con=0
            for j in range(self.matA.matA.dimension):
                con=con+A[i][j]*rVars[j]
            if B[i][0]==-2:
                model.addConstr(con<B[i][1])
            elif B[i][0]==-1:
                model.addConstr(con<=B[i][1])
            elif B[i][0]==0:
                model.addConstr(con==B[i][1])
            elif B[i][0]==1:
                model.addConstr(con>=B[i][1])
            elif B[i][0]==2:
                model.addConstr(con>B[i][1])'''

        model.addConstr(rVars[1]>=0.00001)


        #Produce Objectives

        '''for i in rVars:
            model.setObjective(i)
            model.optimize()
            print(i)
            status = model.Status
            if status==GRB.Status.UNBOUNDED:
                print("UNBOUNDED ",i)
            else:
                if status == GRB.Status.INF_OR_UNBD or \
                   status == GRB.Status.INFEASIBLE  or \
                   status == GRB.Status.UNBOUNDED:
                    print('**The model cannot be solved because it is infeasible or unbounded**')
                else:
                    print("\n\nCounter Example\n\n")
                    for v in model.getVars():
                        print('%s %g' % (v.varName, v.x))
                    print('Obj: %g' % obj.getValue())
                    #return True'''

        model.setObjective(rVars[1])
        model.optimize()
        print(rVars[1])
        status = model.Status
        if status==GRB.Status.UNBOUNDED:
            print("UNBOUNDED ")
        else:
            if status == GRB.Status.INF_OR_UNBD or \
               status == GRB.Status.INFEASIBLE  or \
               status == GRB.Status.UNBOUNDED:
                print('**The model cannot be solved because it is infeasible or unbounded**')
            else:
                print("\n\nCounter Example\n\n")
                for v in model.getVars():
                    print('%s %g' % (v.varName, v.x))
                print('Obj: %g' % obj.getValue())
                return True


            print("------------------------------")

        return False

class Finder:
    def __init__(self,A,B,n1,n2,h,mode):
        self.A=A
        self.B=B
        self.appxM=self.createMatrix(A,B,n1,n2,h,mode)
        self.n=n1+n2

    '''def displayAppxMat(self):
        for i in range(self.n):
            for j in range(self.n):
                print(self.appxM[i][j]," ", end="")
            print('\n')'''

    def createMatrix(self,A,B,n1,n2,h,mode):
        ''' Creates a single matrix based on
        . or +.
        In case of . a roungh approximation is
        done'''

        n=n1+n2
        C=np.zeros((n,n),dtype=np.float128)
        if mode=='+':
            for i in range(n1):
                for j in range(n1):
                    C[i][j]=A[i][j]
            for i in range(n1):
                j2=0
                for j in range(n1,n1+n2):
                    C[i][j]=B[i][j2]
                    j2=j2+1
            for i in range(n1,n1+n2):
                C[i][i]=1
        elif mode=='.':
            I=np.zeros((n1,n1),dtype=np.float128)
            for i in range(n1):
                I[i][i]=1
            A2=h*A
            A2=np.add(I,A2)
            B2=h*B
            for i in range(n1):
                for j in range(n1):
                    C[i][j]=A2[i][j]
            for i in range(n1):
                j2=0
                for j in range(n1,n1+n2):
                    C[i][j]=B2[i][j2]
                    j2=j2+1
            for i in range(n1,n1+n2):
                C[i][i]=1

        return C

    def enumAllBlocks(self):
        '''Enumerates all possible blocks with
        all possible sizes'''
        for m in range(1,self.n+1):
            for k in range(1,self.n+1):
                flt=self.enumBlock(m,k)
                if flt:
                    print("\n\n----",m,k,"----")
                    print(flt)
                    print("====================\n")
                else:
                    print("Finished ",m,",",k)

    def viewRG(self,flt):
        inp=InputFormation(self.appxM,self.n,flt,range(1,len(flt)))
        dyn=Dynamics(inp.formMatA(),self.n,flt)
        rg=RelationGraph([dyn.findConstStr()]+dyn.findVarBlkStr())
        print(rg.isRelationGraph())
        rg.visualizeRelationGraph()

    def enumBlock(self,m,k):
        '''Enumerates possible blocks of
        size m*k'''
        l=m
        b=k
        fltList=[]
        for i in range(self.n-l+1):
            for j in range(self.n-b+1):
                flt=(i,l,j,b)
                inp=InputFormation(self.appxM,self.n,[flt],range(1,len(flt)))
                dyn=Dynamics(inp.formMatA(),self.n,[flt])
                #print(dyn.findConstStr())
                rg=RelationGraph([dyn.findConstStr()]+dyn.findVarBlkStr())
                if (rg.isRelationGraph()):
                    fltList.append(flt)
                    '''print(flt)
                    rg.visualizeRelationGraph()
                    print('\n')'''
        fltList=Finder.filterOutInt(fltList)
        #fltList=[(1,2,3,2),(4,5,6,7),(8,9,10,11),(3,4,5,6)]

        if len(fltList)>0:
            for n in range(len(fltList),0,-1):
                for elem in itertools.combinations(fltList,n):
                    #print(elem)
                    inp=InputFormation(self.appxM,self.n,elem,range(1,len(flt)))
                    dyn=Dynamics(inp.formMatA(),self.n,elem)
                    rg=RelationGraph([dyn.findConstStr()]+dyn.findVarBlkStr())
                    if (rg.isRelationGraph()):
                        return elem


    def filterOutInt(l):
        for elem in l:
            for elem2 in l:
                if elem!=elem2:
                    if (Finder.intersection(elem[0],elem[1],elem2[0],elem2[1]) and Finder.intersection(elem[2],elem[3],elem2[2],elem2[3])):
                        l.remove(elem2)
        #print(l)
        return l

    def intersection(x,lx,y,ly):
        #[[x,lx]] \cap [[y,ly]]
        if (x<=y):
            return (x+lx-1 >= y)
        else:
            return (y+ly-1 >= x)

class BenchmarkSearch:

    def searchTemplate(A,B,n1,n2,h,mode,flt=None):
        fnd=Finder(A,B,n1,n2,h,mode)
        print("-----------")
        if flt==None:
            fnd.enumAllBlocks()
        else:
            fnd.viewRG(flt)

    def madeUp(flt=None):
        BenchmarkSearch.searchTemplate(Benchmarks.Test.A,Benchmarks.Test.B,Benchmarks.Test.n1,Benchmarks.Test.n2,Benchmarks.Test.h,Benchmarks.Test.mode,flt)

    def flightEnvelope(flt=None):
        BenchmarkSearch.searchTemplate(Benchmarks.FlightEnvelope.A,Benchmarks.FlightEnvelope.B,Benchmarks.FlightEnvelope.n1,Benchmarks.FlightEnvelope.n2,Benchmarks.FlightEnvelope.h,Benchmarks.FlightEnvelope.mode,flt)

    def dcConverter(flt=None):
        BenchmarkSearch.searchTemplate(Benchmarks.DCConv.A,Benchmarks.DCConv.B,Benchmarks.DCConv.n1,Benchmarks.DCConv.n2,Benchmarks.DCConv.h,Benchmarks.DCConv.mode,flt)

    def fiveVehiclePlatoon(flt=None):
        BenchmarkSearch.searchTemplate(Benchmarks.FiveVehiclePlatton.A,Benchmarks.FiveVehiclePlatton.B,Benchmarks.FiveVehiclePlatton.n1,Benchmarks.FiveVehiclePlatton.n2,Benchmarks.FiveVehiclePlatton.h,Benchmarks.FiveVehiclePlatton.mode,flt)

    def tenVehiclePlatoon(flt=None):
        BenchmarkSearch.searchTemplate(Benchmarks.TenVehiclePlatton.A,Benchmarks.TenVehiclePlatton.B,Benchmarks.TenVehiclePlatton.n1,Benchmarks.TenVehiclePlatton.n2,Benchmarks.TenVehiclePlatton.h,Benchmarks.TenVehiclePlatton.mode,flt)

    def coOPVehicles(flt=None):
        BenchmarkSearch.searchTemplate(Benchmarks.CoOPVehicles.A,Benchmarks.CoOPVehicles.B,Benchmarks.CoOPVehicles.n1,Benchmarks.CoOPVehicles.n2,Benchmarks.CoOPVehicles.h,Benchmarks.CoOPVehicles.mode,flt)

    def pkpd(flt=None):
        BenchmarkSearch.searchTemplate(Benchmarks.PKPD.A,Benchmarks.PKPD.B,Benchmarks.PKPD.n1,Benchmarks.PKPD.n2,Benchmarks.PKPD.h,Benchmarks.PKPD.mode,flt)

    def pkpd2(flt=None):
        BenchmarkSearch.searchTemplate(Benchmarks.PKPD2.A,Benchmarks.PKPD2.B,Benchmarks.PKPD2.n1,Benchmarks.PKPD2.n2,Benchmarks.PKPD2.h,Benchmarks.PKPD2.mode,flt)

    def spaceCraftRndzvs(flt=None):
        BenchmarkSearch.searchTemplate(Benchmarks.SpaceCraftRndzvs.A,Benchmarks.SpaceCraftRndzvs.B,Benchmarks.SpaceCraftRndzvs.n1,Benchmarks.SpaceCraftRndzvs.n2,Benchmarks.SpaceCraftRndzvs.h,Benchmarks.SpaceCraftRndzvs.mode,flt)

    def holesCXc(flt=None):
        BenchmarkSearch.searchTemplate(Benchmarks.HolesCXc.A,Benchmarks.HolesCXc.B,Benchmarks.HolesCXc.n1,Benchmarks.HolesCXc.n2,Benchmarks.HolesCXc.h,Benchmarks.HolesCXc.mode,flt)

    def holesPDp(flt=None):
        BenchmarkSearch.searchTemplate(Benchmarks.HolesPDp.A,Benchmarks.HolesPDp.B,Benchmarks.HolesPDp.n1,Benchmarks.HolesPDp.n2,Benchmarks.HolesPDp.h,Benchmarks.HolesPDp.mode,flt)

    def holesPXp(flt=None):
        BenchmarkSearch.searchTemplate(Benchmarks.HolesPXp.A,Benchmarks.HolesPXp.B,Benchmarks.HolesPXp.n1,Benchmarks.HolesPXp.n2,Benchmarks.HolesPXp.h,Benchmarks.HolesPXp.mode,flt)

    def motorTransmission1(flt=None):
        BenchmarkSearch.searchTemplate(Benchmarks.MotorTransmission1.A,Benchmarks.MotorTransmission1.B,Benchmarks.MotorTransmission1.n1,Benchmarks.MotorTransmission1.n2,Benchmarks.MotorTransmission1.h,Benchmarks.MotorTransmission1.mode,flt)

    def motorTransmission2(flt=None):
        BenchmarkSearch.searchTemplate(Benchmarks.MotorTransmission2.A,Benchmarks.MotorTransmission2.B,Benchmarks.MotorTransmission2.n1,Benchmarks.MotorTransmission2.n2,Benchmarks.MotorTransmission2.h,Benchmarks.MotorTransmission2.mode,flt)

class Benchmarks:

    def createMatrix(A,B,n1,n2,h,mode):
        ''' Creates a single matrix based on
        . or +.
        In case of . a roungh approximation is
        done'''

        n=n1+n2
        C=np.zeros((n,n),dtype=np.float128)
        if mode=='+':
            for i in range(n1):
                for j in range(n1):
                    C[i][j]=A[i][j]
            for i in range(n1):
                j2=0
                for j in range(n1,n1+n2):
                    C[i][j]=B[i][j2]
                    j2=j2+1
            for i in range(n1,n1+n2):
                C[i][i]=1
        elif mode=='.':
            I=np.zeros((n1,n1),dtype=np.float128)
            for i in range(n1):
                I[i][i]=1
            A2=h*A
            A2=np.add(I,A2)
            B2=h*B
            for i in range(n1):
                for j in range(n1):
                    C[i][j]=A2[i][j]
            for i in range(n1):
                j2=0
                for j in range(n1,n1+n2):
                    C[i][j]=B2[i][j2]
                    j2=j2+1
            for i in range(n1,n1+n2):
                C[i][i]=1

        return C

    class Test:
        '''A=np.array([
        [1,0,0,1],
        [0,0,0,0],
        [0,0,1,0],
        [0,0,0,1],
        ])
        B=np.array([
        [1,0],
        [0,1],
        [0,0],
        [0,3],
        ])'''

        A=np.array([
        [1,0,0,0,0,0],
        [1,0,0,0,0,0],
        [1,0,0,0,0,0],
        [1,0,0,0,0,0],
        [1,0,0,0,0,0],
        [1,0,0,0,0,0],
        ])
        B=np.array([
        ])
        n1=6
        n2=0
        h=0.01
        mode='.'

    class FlightEnvelope:
        A=np.array([
        [0,0,0,1,0,0,0,0,0,0,0,0],
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
        ])
        B=np.array([
        [0,0,0,0],
        [0,0,0,0],
        [0,0,0,0],
        [0,0,0,0],
        [0,0,0,0],
        [1,0,0,0],
        [0,0,0,0],
        [0,0,0,0],
        [0,0,0,0],
        [0,1,0,0],
        [0,0,1,0],
        [0,0,0,1],
        ])
        n1=12
        n2=4
        h=0.01
        mode='.'

    class DCConv:
        #IMPORTANT: Values are inaccuarate
        A=np.array([
        [0,0],
        [0,-5],
        ])
        B=np.array([
        [5],
        [0]
        ])
        n1=2
        n2=1
        h=0.01
        mode='.'

    class FiveVehiclePlatton:
        A=np.array([
        [0,1,0,0,0,0,0,0,0,0,0,0,0,0,0],
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
        B=np.array([
        ])
        n1=15
        n2=0
        h=0.01
        mode='.'

    class TenVehiclePlatton:
        A=np.array([
        [0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
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
        ])
        B=np.array([
        ])
        n1=30
        n2=0
        h=0.01
        mode='.'

    class CoOPVehicles:
        A=np.array([
        [0,1,0,0,0,0,0,0,0],
        [0,0,-1,0,0,0,0,0,0],
        [1.6050,4.8680,-3.5754,-0.8198,0.4270,-0.0450,-0.1942,0.3626,-0.0946],
        [0,0,0,0,1,0,0,0,0],
        [0,0,1,0,0,-1,0,0,0],
        [0.8718,3.8140,-0.0754,1.1936,3.6258,-3.2396,-0.5950,0.1294,-0.0796],
        [0,0,0,0,0,0,0,1,0],
        [0,0,0,0,0,1,0,0,-1],
        [0.7132,3.5730,-0.0964,0.8472,3.2568,-0.0876,1.2726,3.0720,-3.1356]
        ])
        B=np.array([
        [0],
        [1],
        [0],
        [0],
        [0],
        [0],
        [0],
        [0],
        [0],
        ])
        n1=9
        n2=1
        h=0.01
        mode='.'
        '''A=np.array([
        [0,1,0,0,0,0,0,0,0],
        [0,0,-1,0,0,0,0,0,0],
        [1.6050,4.8680,-3.5754,0,0,0,0,0,0],
        [0,0,0,0,1,0,0,0,0],
        [0,0,1,0,0,-1,0,0,0],
        [0,0,0,1.1936,3.6258,-3.2396,0,0,0],
        [0,0,0,0,0,0,0,1,0],
        [0,0,0,0,0,1,0,0,-1],
        [0.7132,3.5730,-0.0964,0.8472,3.2568,-0.0876,1.2726,3.0720,-3.1356],
        ])
        B=np.array([
        [0],
        [1],
        [0],
        [0],
        [0],
        [0],
        [0],
        [0],
        [0],
        ])'''

    class PKPD:
        #With weight of 25 Kg
        weight=25
        v1=458.4*weight
        k10=0.1527*pow(weight,-0.3)
        k12=0.114
        k13=0.0419
        k21=0.055
        k31=0.0033
        A=np.array([
        [-(k10+k12+k13),k12,k13],
        [k21,-k21,0],
        [k31,0,-k31]
        ])
        B=np.array([
        [1/v1],
        [0],
        [0]
        ])
        n1=3
        n2=1
        h=0.01
        mode='.'

    class PKPD2:
        #With weight of 25 Kg
        weight=25
        v1=458.4*weight
        k10=0.1527*pow(weight,-0.3)
        k12=0.114
        k13=0.0419
        k21=0.055
        k31=0.0033
        kd=40 #with Td=20s
        A=np.array([
        [-(k10+k12+k13),k12,k13,0],
        [k21,-k21,0,0],
        [k31,0,-k31,0],
        [kd,0,0,-kd]
        ])
        B=np.array([
        [1/v1],
        [0],
        [0],
        [0]
        ])
        n1=4
        n2=1
        h=0.01
        mode='.'

    class SpaceCraftRndzvs:
        #With weight of 25 Kg
        mu=3.986e14
        r=42164
        mc=500
        n=math.sqrt(mu/(pow(r,3)))
        A=np.array([
        [0,0,1,0],
        [0,0,0,1],
        [3*n*n,0,0,2*n],
        [0,0,-2*n,0]
        ])
        B=np.array([
        [0,0],
        [0,0],
        [1/mc,0],
        [0,1/mc]
        ])
        n1=4
        n2=2
        h=0.01
        mode='.'

    class HolesCXc:
        A=np.array([
        [0,1,0,-1,0],
        [0,0,0.996,0,0],
        [0,0,-0.9999999525,0,0],
        [0,0,0,0,0],
        [0,0,0,49.99999763,-49.99999763]
        ])
        B=np.array([
        [0,-0.00872664619,-0.00872664619,0,0],
        [0,0,0,0.0003,0],
        [0,0,0,0.07499999644,0],
        [0.999999992,0,0,0,0],
        [0,0.4363322888,0.4363322888,0,0]
        ])
        n1=5
        n2=5
        h=0.01
        mode='.'

    class HolesPDp:
        A=np.array([
        [0.08087,0],
        [0,0.08087]
        ])
        B=np.array([
        ])
        n1=2
        n2=0
        h=0.01
        mode='.'

    class HolesPXp:
        A=np.array([
        [-162.1272879,0,162.1272879,-409.7154124,0,0],
        [1,0,0,0,0,0],
        [78.14959139,0,-78.14959139,269.6273393,0,0],
        [0,0,1,0,0,0],
        [0,0,0,0,-92.4134822,0],
        [0,0,0,0,1,0],
        ])
        B=np.array([
        [157.5738419,157.5738419],
        [0,0],
        [-75.95760351,-75.95760351],
        [0,0],
        [-51.3265211,51.3265211],
        [0,0]
        ])
        n1=6
        n2=2
        h=0.01
        mode='.'

    class MotorTransmission1:
        '''Needs change of matrices'''

        def createMatrix(A,B,n1,n2,h,mode):
            ''' Creates a single matrix based on
            . or +.
            In case of . a roungh approximation is
            done'''

            n=n1+n2
            C=np.zeros((n,n),dtype=np.float128)
            if mode=='+':
                for i in range(n1):
                    for j in range(n1):
                        C[i][j]=A[i][j]
                for i in range(n1):
                    j2=0
                    for j in range(n1,n1+n2):
                        C[i][j]=B[i][j2]
                        j2=j2+1
                for i in range(n1,n1+n2):
                    C[i][i]=1
            elif mode=='.':
                I=np.zeros((n1,n1),dtype=np.float128)
                for i in range(n1):
                    I[i][i]=1
                A2=h*A
                A2=np.add(I,A2)
                B2=h*B
                for i in range(n1):
                    for j in range(n1):
                        C[i][j]=A2[i][j]
                for i in range(n1):
                    j2=0
                    for j in range(n1,n1+n2):
                        C[i][j]=B2[i][j2]
                        j2=j2+1
                for i in range(n1,n1+n2):
                    C[i][i]=1

            return C

        ms=3.2
        Rs=0.08
        jg2=0.09
        Tf=1
        A=np.array([
        [0,0,0,0,0],
        [0,0,0,0,0],
        [1,0,0,0,0],
        [0,1,0,0,0],
        [0,0,0,0,0],
        ])
        B=np.array([
        [1/ms,0],
        [0,-(Rs*Tf)/jg2],
        [0,0],
        [0,0],
        [0,0]
        ])
        '''fi=np.array([
        [0,-(Rs*Tf)/jg2,0,0,0]
        ])'''
        n1=5
        n2=2
        h=0.01
        mode='.'

    class MotorTransmission2:
        '''Needs change of matrices'''

        ms=3.2
        A=np.array([
        [0,0,0,0,0],
        [0,0,0,0,0],
        [0,0,1,0,0],
        [0,1,0,1,0],
        [ms,ms,0,0,1],
        ])
        B=np.array([])
        n1=5
        n2=0
        h=0.01
        mode='.'

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


    def experimentalGenerator5Help(par):
        A=par[0][0]
        n=par[0][1]
        r=par[1][0]
        c=par[1][1]
        flt=[(r,c)]
        #if A[r][c]!=0:
        inp=InputFormation(A,n,flt,range(1,len(flt)),range(1,len(flt)),5)
        dynA=Dynamics(inp.formMatA(),n,flt)
        rg=RelationGraph([dynA.findConstStr()]+dynA.findVarBlkStr())
        if (rg.isRelationGraph()):
            return (r,c)

    def createMatrix(A,B,n1,n2,h,mode):
        n=n1+n2
        C=np.zeros((n,n),dtype=np.float128)
        if mode=='+':
            for i in range(n1):
                for j in range(n1):
                    C[i][j]=A[i][j]
            for i in range(n1):
                j2=0
                for j in range(n1,n1+n2):
                    C[i][j]=B[i][j2]
                    j2=j2+1
            for i in range(n1,n1+n2):
                C[i][i]=1
        elif mode=='.':
            I=np.zeros((n1,n1),dtype=np.float128)
            for i in range(n1):
                I[i][i]=1
            A2=h*A
            A2=np.add(I,A2)
            B2=h*B
            for i in range(n1):
                for j in range(n1):
                    C[i][j]=A2[i][j]
            for i in range(n1):
                j2=0
                for j in range(n1,n1+n2):
                    C[i][j]=B2[i][j2]
                    j2=j2+1
            for i in range(n1,n1+n2):
                C[i][i]=1

        return C


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
        A1=np.array([
        [0,0,0,1,0,0,0,0,0,0,0,0],
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
        ])
        B1=np.array([
        [0,0,0,0],
        [0,0,0,0],
        [0,0,0,0],
        [0,0,0,0],
        [0,0,0,0],
        [1,0,0,0],
        [0,0,0,0],
        [0,0,0,0],
        [0,0,0,0],
        [0,1,0,0],
        [0,0,1,0],
        [0,0,0,1],
        ])
        n1=12
        n2=4
        h=1
        A=Examples.createMatrix(A1,B1,n1,n2,h,'.')
        n=n1+n2
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


    def flightEnvelope():
        A=np.array([
        [0,0,0,1,0,0,0,0,0,0,0,0],
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
        ])
        B=np.array([
        [0,0,0,0],
        [0,0,0,0],
        [0,0,0,0],
        [0,0,0,0],
        [0,0,0,0],
        [1,0,0,0],
        [0,0,0,0],
        [0,0,0,0],
        [0,0,0,0],
        [0,1,0,0],
        [0,0,1,0],
        [0,0,0,1],
        ])
        n1=12
        n2=4
        f=[(0,12),(0,14),(0,15),(2,13),(8,12),(8,15)]
        v=[5,5,6,6,7,8] #Yet to be decided
        T=20 #Yet to be figured out
        h=1 #Yet to be figured out
        Ca=np.array([
        [1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
        ])
        Cb=np.array([
        [-1, 2],
        ])
        l=1
        initset=[(3,4),(3,4),(3,5),(3,4),(3,4),(3,5),(-0.50,0.50),(-0.50,0.50),(-0.50,0.50),(-1,1),(-1,1),(-1,1),(5,10),(-0.1,0.1),(-0.1,0.1),(-0.1,0.1)] #Yet to figured out, says uncertain
        #initset=[(0,5),(-0.08,0.08),(0,5),(-0.01,0.01),(0,30)]
        Examples.exampleTemplate(A,B,'.',h,n1,n2,f,v,T,initset,"Flight Envelope",Ca,Cb,l)

    def holes():
        A=np.array([[1,0.00400000019,0,-0.00400000019,0],
                    [0,1,0.003984000189,0,0],
                    [0,0,0.996,0,0],
                    [0,0,0,1,0],
                    [0,0,0,0.2,0.8],
                    ])
        B=np.array([[0,-0.00003490658642,-0.00003490658642,0,0],
                    [0,0,0,0.000001200000057,0],
                    [0,0,0,0.0003,0],
                    [0.004000000158,0,0,0,0],
                    [0,0.001745329238,0.001745329238,0,0],
                    ])
        n1=5
        n2=5
        f=[(0,6),(0,7),(4,6),(4,7)]
        v=[5,5,6,7]
        T=7
        h=1
        Ca=np.array([
        [1,0,0,0,0,0,0,0,0,0],
        [0,1,0,0,0,0,0,0,0,0],
        [0,0,1,0,0,0,0,0,0,0],
        [0,0,0,1,0,0,0,0,0,0],
        [0,0,0,0,1,0,0,0,0,0],
        [0,0,0,0,0,1,0,0,0,0],
        [0,0,0,0,0,0,1,0,0,0],
        [0,0,0,0,0,0,0,1,0,0],
        [0,0,0,0,0,0,0,0,1,0],
        [0,0,0,0,0,0,0,0,0,1],
        ])
        Cb=np.array([
        [1,1.39],
        [1,1.39],
        [1,1.39],
        [1,1.39],
        [1,1.39],
        [1,1.39],
        [1,1.39],
        [1,1.39],
        [1,1.39],
        [1,1.39],
        ])
        l=10
        initset=[(0,0),(0,0),(0,0),(0,0),(0,0),(-100,100),(0,0),(0,0),(-100,100),(-100,100)]
        Examples.exampleTemplate(A,B,'+',h,n1,n2,f,v,T,initset,"Holes Controller - XC",Ca,Cb,l)

    def transformer():
        A=np.array([[0,0,0],
                    [0,0,0],
                    [0,0,0.00003],
                    ]) #Yet to figured out
        B=np.array([[12],
                    [0],
                    [0]])
        n1=3
        n2=1
        f=[(0,3)]
        #v=[50] #Yet to be decided
        v=[50]
        T=20 #Yet to be figured out
        h=0.01 #Yet to be figured out
        Ca=np.array([[0,1,0,0],
        ])
        Cb=np.array([[1,0.0000001],
        ])
        '''Ca=np.array([[1,0,0,0],
        ])
        Cb=np.array([[1,10],
        ])'''
        l=1
        initset=[(0,0.4),(-5,0),(0,0.4),(1,1)] #Yet to figured out
        Examples.exampleTemplate(A,B,'.',h,n1,n2,f,v,T,initset,"Transformer Isolated DC-DC Converter",Ca,Cb,l)

    def exampleTemplate(At,Bt,mode,h,n1,n2,f,v,T,initset,strName,Ca,Cb,l):

        A=Examples.createMatrix(At,Bt,n1,n2,h,mode)
        n=n1+n2
        print("----Final Matrix----")
        for i in range(n):
            for j in range(n):
                print(A[i][j], end="  ")
            print("\n")
        print("---------------------")
        gurobi_time=0
        algo_time=0
        start_time=time.time()
        inpA=InputFormation(A,n,f,v,initset,T)
        dynA=Dynamics(inpA.formMatA(),n,f)
        rg=RelationGraph([dynA.findConstStr()]+dynA.findVarBlkStr())
        if (rg.isRelationGraph()==True):
            for i in range(1,int(T/h)+1):
                print("\n\n\n\n==============",i)
                dynRepA=DynamicsRepOp(dynA)
                reach=CalcReachableSet(dynRepA,i,initset,f,inpA.perturbations())
                reach.matAk.displayA()
                flag=reach.isFeasible(Ca,Cb,l)
                end_time=time.time() - start_time
                if flag==True:
                    dynA.displayA()
                    print("UNSAFE!!")
                    print(end_time)
                    exit()
        print("\n\nSafe System!\n\n")
        print(end_time)



BenchmarkSearch.flightEnvelope()
