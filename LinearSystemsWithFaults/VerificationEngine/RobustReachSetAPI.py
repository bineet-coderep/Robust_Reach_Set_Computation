'''
Author: Bineet Ghosh, under supervision of Prof. Sridhar
Email: ghosh.bineet22@gmail.com

Code Freeze Date: 

Based on the paper:
'Robust Reachable Set: Accounting for Uncertainties in Linear Dynamical Systems'
by Bineet Ghosh and Parasara Sridhar Duggirala

 - Given a Linear Dynamical System, this engine can find all places where fault can
be induced such that the verification can be performed efficiently.

- Given a Linear Dynamical System with uncertainties (that can either be constant
with time or varying), this engine can verfify the safety of such systems efficiently

Please refer to the README.md to use or modify these APIs.

The main Verification Engine is: uncertain-system-verifier.py
which uses these APIs

Link to the repository: <>
'''

import numpy as np
import itertools
import sys
import time
import math
from gurobipy import *
from numpy import linalg as LA
import xml.etree.ElementTree as ET

class Structure:
    '''Represents Support of a matrix
    '''

    def __init__(self,l,n):
        self.omegas=l #List tuple of indices: (startrow,length,startcol,breadth)
        self.size=n #The support is of dimension n*n
        '''creates an actual matrix Representing this structure
        1 represents the set R, and 0 otherwise '''
        self.matrix=np.zeros((self.size,self.size),dtype=int)
        for omega in self.omegas:
            for i in range(omega[0],omega[0]+omega[1]):
                for j in range(omega[2],omega[2]+omega[3]):
                    self.matrix[i][j]=1 #1 Represents the set R

    def visualizeStr(self):
        '''displays the support in form a Matrix
        '''
        for i in range(0,self.size):
            for j in range(0,self.size):
                print (self.matrix[i][j], " ",end="")
            print ("\n")

class BlockStructure(Structure):
    '''Represents Block Structured Matrices
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
        #Check if the product of this support and the given support is zero.
        return (BlockStructure.intersection(self.getColStart(),self.getColBreadth(),mat.getRowStart(),mat.getRowLength())==False and BlockStructure.intersection(mat.getColStart(),mat.getColBreadth(),self.getRowStart(),self.getRowLength())==False)

    def getFaultyCells(self):
        #Returns the list of cells contained in the block
        l=[]
        for i in range(self.getRowStart(),self.getRowStart()+self.getRowLength()):
            for j in range(self.getColStart(),self.getColStart()+self.getColBreadth()):
                l.append((i,j))
        return l


    #def generateAll()

class StructLME:
    '''Represents the list of Structures
    corresponding to all the variables i.e. N_i
    and the constant structure N_0 too.
    '''
    def __init__(self, l):
        self.vertices=l #l[0]=N_0 and rest l[i]=N_i
        self.numOfVars=len(l)-1

    def isCondSAT(self):
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
        '''Checks if support of N_0 satisfies the conditions mentioned in the paper
        including closure condition'''
        U=np.ones((self.vertices[1].size,self.vertices[1].size),dtype=int)

        for elem in self.vertices[1:]:
            U=StructLME.intersectionOfBlocks(U,StructLME.findStrN0ForOne(elem),self.vertices[1].size)

        for i in range(0,self.vertices[1].size):
            for j in range(0,self.vertices[1].size):
                if (U[i][j]==0 and self.vertices[0].matrix[i][j]!=0):
                    return False

        '''print(U)
        print(StructLME.isStrClosed(U,self.vertices[1].size))
        exit(0)'''
        return StructLME.isStrClosed(U,self.vertices[1].size)
        #return StructLME.isStrClosed(self.vertices[0].matrix,self.vertices[0].size)

    @staticmethod
    def isStrClosed(U,n):
        U2=np.matmul(U,U)
        for i in range(0,n):
            for j in range(0,n):
                if U[i][j]==0 and U2[i][j]!=0:
                    return False

        #print(U)
        return True

    def visualizeStructLME(self):
        #Displays the list of Supports
        self.vertices[0].visualizeStr() #Structure of N_0
        print()
        for elem in self.vertices[1:]:
            elem.visualizeStr() #BlockStructure of N_i
            print()

    @staticmethod
    def findStrN0ForOne(blk):
        '''Returns a super support for N_0 such that blk
        satisfies the conditions mentioned in the paper
        '''
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

        return (StructLME.intersectionOfBlocks(O,M,blk.size))

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
        '''Finds out the Support of N_0
        of the given dynamics.
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
        '''Finds out the N_0
        of the given dynamics.
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
            for i in range(v.getRowStart(),v.getRowStart()+v.getRowLength()):
                for j in range(v.getColStart(),v.getColStart()+v.getColBreadth()):
                    ar[i][j]=self.inpA[i][j]
            l.append(ar)
        return l

class DynamicsRepOp:
    '''Represents the LME of the given dynamics
    '''

    def __init__(self,N0,varNs,s):
        self.N0=N0
        self.varNs=varNs
        self.dimension=s

    def displayA(self):
        for i in range(0,self.matA.dimension):
            for j in range(0,self.matA.dimension):
                print (N0[i][j],"  ",end="")
            print("\n")

        print("\n")
        for var in varNs:
            for i in range(0,self.dimension):
                for j in range(0,self.dimension):
                    print (var[i][j],"  ",end="")
                print("\n")
            print("\n")

    def power(self, k):
        '''Calculates the power k
        '''
        M0=[]
        varMs=[]
        M0[:]=self.N0[:]
        varMs[:]=self.varNs[:]
        for x in range(0,k-1):
            for i in range(0,len(varMs)):
                varMs[i]=np.add(np.matmul(self.varNs[i],M0),np.matmul(self.N0,varMs[i]))
            M0=np.matmul(self.N0,M0)

        return [M0]+varMs

    def powerTimeVarying(self, k, nNoises):
        '''Calculates the power k,
        considering the non constant
        noises
        '''

        M0=[]
        varMs=[]
        M0[:]=self.N0[:]
        for i in range(0,len(self.varNs)):
            varMs.append([self.varNs[i]])

        '''for i in range(0,len(self.varNs)):
            print (varMs[i][0])'''

        crntTime=1
        for pow in range(0,k):
            #Multiply two LMEs
            if pow%nNoises!=0 or pow==0:
                for i in range(0,len(self.varNs)):
                    for j in range(0,len(varMs[i])):
                        varMs[i][j]=np.matmul(varMs[i][j],self.N0)
                for i in range(0,len(self.varNs)):
                    varMs[i][crntTime-1]=np.add(varMs[i][crntTime-1],np.matmul(M0,self.varNs[i]))
            else:
                crntTime=crntTime+1
                for i in range(0,len(self.varNs)):
                    for j in range(0,len(varMs[i])):
                        varMs[i][j]=np.matmul(varMs[i][j],self.N0)
                for i in range(0,len(self.varNs)):
                    varMs[i].append(np.matmul(M0,self.varNs[i]))
            M0=np.matmul(M0,self.N0)

            '''for i in range(0,len(self.varNs)):
                print (len(varMs[i]))'''

        return [M0]+varMs

class InputFormation:

    def __init__(self,A, n, v):
        self.dynA=A
        self.dimension=n
        self.vars=v
        #self.variance=vrnc

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

    '''def perturbations(self):
        pr=[]
        for i in range(len(self.vars)):
            mn=self.minimumInCells(i)
            mx=self.maximumInCells(i)
            t1=mn-((mn*self.variance[i])/100)
            t2=mx+((mx*self.variance[i])/100)
            print((min(t1,t2),max(t1,t2)))
            pr.append((min(t1,t2),max(t1,t2)))
        #pr=[(-0.0000001,0)]
        return pr

    def minimumInCells(self,v):
        e=self.vars[v]
        mn=9999999
        for i in range(e[0],e[0]+e[1]):
            for j in range(e[2],e[2]+e[3]):
                if self.dynA[i][j]<mn:
                    mn=self.dynA[i][j]
        return mn

    def maximumInCells(self,v):
        e=self.vars[v]
        mx=-9999999
        for i in range(e[0],e[0]+e[1]):
            for j in range(e[2],e[2]+e[3]):
                if self.dynA[i][j]>mx:
                    mx=self.dynA[i][j]
        return mx'''

class CalcReachableSet:

    def __init__(self,A,k,initset,flt,prtb,unsafeV,unsafeP):
        self.matA=A
        self.matAk=A.power(k)
        self.initialSet=initset
        self.faults=flt
        self.perturbations=prtb
        self.unsafeVariable=unsafeV
        self.unsafeParam=unsafeP

    def getAllFaultyCells(self):
        l=[]
        for e in self.faults:
            for i in range(e[0],e[0]+e[1]):
                for j in range(e[2],e[2]+e[3]):
                    l.append((i,j))
        return l

    def indexOfBlk(self, i,j):
        for ind in range(len(self.faults)):
            e=self.faults[ind]
            if i>=e[0] and i<e[0]+e[1]:
                if j>=e[2] and j<e[2]+e[3]:
                    return ind

    def isFeasible(self):

        flag=True
        model = Model("qp")
        #print(self.perturbations[0])
        #exit(0)
        #self.perturbations[0]=(0,0.6)
        #print(self.perturbations[0])
        #exit(0)
        #model.Params.BarHomogeneous=1
        #model.Params.NumericFocus=3

        #create variables
        faultVars=[]
        for i in self.faults:
            name="fault "+str(i)
            faultVars.append(model.addVar(-GRB.INFINITY,GRB.INFINITY,name=name,vtype='C'))
        dynVars=[]
        for i in range(self.matA.dimension):
            name="alpha"+str(i)
            dynVars.append(model.addVar(-GRB.INFINITY,GRB.INFINITY,name=name,vtype='C'))

        #Add initial Set and Perturbations constraints
        for i in range(self.matA.dimension):
            model.optimize()
            name="Init-C"+str(i)
            model.addConstr(dynVars[i]>=self.initialSet[i][0],name+".1")
            model.addConstr(dynVars[i]<=self.initialSet[i][1],name+".2")

        no=0
        for i in faultVars:
            name="Fault-C"+str(no)
            model.addConstr(i>=self.perturbations[no][0],name+".1")
            model.addConstr(i<=self.perturbations[no][1],name+".2")
            no+=1

        #Generate Final Equations
        rVars=[]
        for i in range(self.matA.dimension):
            obj=0
            for j in range(self.matA.dimension):
                if (i,j) in self.getAllFaultyCells():
                    indx=self.indexOfBlk(i,j)
                    obj=obj+((self.matAk[0][i][j]+self.matAk[indx+1][i][j]*faultVars[indx])*dynVars[j])
                    #obj=obj+(self.matAk.inpA[i][j]*faultVars[self.indexOfBlk(i,j)]*dynVars[j])
                else:
                    obj=obj+(self.matAk[0][i][j]*dynVars[j])
            rVars.append(obj)


        model.addConstr(rVars[self.unsafeVariable]>=self.unsafeParam)

        #Produce Objectives

        model.setObjective(rVars[self.unsafeVariable])
        model.optimize()
        #print(rVars[0])
        #exit()
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

class CalcReachableSetTimeVarying:

    def __init__(self,A,k,initset,flt,prtb,unsafeV,unsafeP,t):
        self.matA=A
        self.matAk=A.powerTimeVarying(k,t)
        self.initialSet=initset
        self.faults=flt
        self.perturbations=prtb
        self.unsafeVariable=unsafeV
        self.unsafeParam=unsafeP

    def getAllFaultyCells(self):
        l=[]
        for e in self.faults:
            for i in range(e[0],e[0]+e[1]):
                for j in range(e[2],e[2]+e[3]):
                    l.append((i,j))
        return l

    def indexOfBlk(self, i,j):
        for ind in range(len(self.faults)):
            e=self.faults[ind]
            if i>=e[0] and i<e[0]+e[1]:
                if j>=e[2] and j<e[2]+e[3]:
                    return ind

    def isFeasible(self):

        flag=True
        model = Model("qp")
        #print(self.perturbations[0])
        #exit(0)
        #self.perturbations[0]=(0,0.6)
        #print(self.perturbations[0])
        #exit(0)
        #model.Params.BarHomogeneous=1
        #model.Params.NumericFocus=3

        #create variables
        #print(self.matAk[1])
        #exit(0)
        faultVars=[]
        timeV=len(self.matAk[1])
        for i in self.faults:
            l=[]
            for j in range(timeV):
                name="fault "+str(i)+"-"+str(j)
                l.append(model.addVar(-GRB.INFINITY,GRB.INFINITY,name=name,vtype='C'))
            faultVars.append(l)
        dynVars=[]
        for i in range(self.matA.dimension):
            name="alpha"+str(i)
            dynVars.append(model.addVar(-GRB.INFINITY,GRB.INFINITY,name=name,vtype='C'))

        #Add initial Set and Perturbations constraints
        for i in range(self.matA.dimension):
            model.optimize()
            name="Init-C"+str(i)
            model.addConstr(dynVars[i]>=self.initialSet[i][0],name+".1")
            model.addConstr(dynVars[i]<=self.initialSet[i][1],name+".2")

        no=0
        for i in faultVars:
            no2=0
            for j in i:
                name="FltC"+str(no)+"-"+str(j)
                model.addConstr(j>=self.perturbations[no][no2][0],name+".1")
                model.addConstr(j<=self.perturbations[no][no2][1],name+".2")
                no2=no2+1
            no+=1

        #Generate Final Equations
        rVars=[]
        '''for i in range(self.matA.dimension):
            obj=0
            for j in range(self.matA.dimension):
                if (i,j) in self.getAllFaultyCells():
                    indx=self.indexOfBlk(i,j)
                    obj=obj+((self.matAk[0][i][j]+self.matAk[indx+1][i][j]*faultVars[indx])*dynVars[j])
                    #obj=obj+(self.matAk.inpA[i][j]*faultVars[self.indexOfBlk(i,j)]*dynVars[j])
                else:
                    obj=obj+(self.matAk[0][i][j]*dynVars[j])
            rVars.append(obj)'''

        for i in range(self.matA.dimension):
            obj=0
            for j in range(self.matA.dimension):
                if (i,j) in self.getAllFaultyCells():
                    indx=self.indexOfBlk(i,j)
                    #obj=obj+((self.matAk[0][i][j]+self.matAk[indx+1][i][j]*faultVars[indx])*dynVars[j])
                    obj=obj+self.matAk[0][i][j]
                    for k in range(timeV):
                        obj=obj+(self.matAk[indx+1][k][i][j]*faultVars[indx][k]*dynVars[j])
                else:
                    obj=obj+(self.matAk[0][i][j]*dynVars[j])
            rVars.append(obj)


        model.addConstr(rVars[self.unsafeVariable]<=self.unsafeParam)

        #Produce Objectives

        model.setObjective(rVars[self.unsafeVariable])
        model.optimize()
        '''print("FS")
        print(rVars[2])
        exit()'''
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

class CalcReachableSetFaultFree:

    def __init__(self,A,n,k,initset,unsafeV,unsafeP):
        self.matA=A
        self.dimension=n
        self.matAk=LA.matrix_power(A,k)
        self.initialSet=initset
        self.unsafeVariable=unsafeV
        self.unsafeParam=unsafeP

    def isFeasible(self):

        flag=True
        model = Model("lp")
        #print(self.perturbations[0])
        #exit(0)
        #self.perturbations[0]=(0,0.6)
        #print(self.perturbations[0])
        #exit(0)
        #model.Params.BarHomogeneous=1
        #model.Params.NumericFocus=3

        #create variables
        dynVars=[]
        for i in range(self.dimension):
            name="alpha"+str(i)
            dynVars.append(model.addVar(-GRB.INFINITY,GRB.INFINITY,name=name,vtype='C'))

        #Add initial Set and Perturbations constraints
        for i in range(self.dimension):
            model.optimize()
            name="Init-C"+str(i)
            model.addConstr(dynVars[i]>=self.initialSet[i][0],name+".1")
            model.addConstr(dynVars[i]<=self.initialSet[i][1],name+".2")

        #Generate Final Equations
        rVars=[]
        for i in range(self.dimension):
            obj=0
            for j in range(self.dimension):
                obj=obj+(self.matAk[i][j]*dynVars[j])
            rVars.append(obj)

        model.addConstr(rVars[self.unsafeVariable]<=self.unsafeParam)

        #Produce Objectives

        model.setObjective(rVars[self.unsafeVariable])
        model.optimize()
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

class Benchmarks:

    def createMatrix(A,B,h,mode):
        ''' Creates a single matrix based on
        . or +.
        In case of . a roungh approximation is
        done'''

        n1=np.size(A,0)
        if (np.size(B)>0):
            n2=np.size(B,1)
        else:
            n2=0
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

    class IllustExample:
        A=np.array([
        [3,0,0,0,0,2,4],
        [1,2,3,0,0,2.9,2.9],
        [8,1,2,0,0,2.9,2.9],
        [7,0,0,8,2,3.9,3.9],
        [8,0,0,3,7,3.9,3.9],
        [0,0,0,0,0,6,3],
        [0,0,0,0,0,2,1],
        ])
        B=np.array([
        [3],
        [0],
        [0],
        [9],
        [11],
        [0],
        [0],
        ])
        h=0.01
        mode='+'

    class IllustExample2:
        A=np.array([
        [3,2.9,3.9],
        [0,7,0],
        [0,0,2]
        ])
        B=np.array([
        [2],
        [0],
        [0],
        ])
        h=0.01
        mode='+'

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
        ])

        '''
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
        h=0.01
        mode='.'

    class CoOPVehiclesI:
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
        h=0.01
        mode='.'

    class CoOPVehiclesII:
        A=np.array([
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
        ])
        h=0.01
        mode='.'

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
        h=0.01
        mode='.'

    class SpaceCraftRndzvs:
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
        h=0.01
        mode='.'

    class HolesPDp:
        A=np.array([
        [0.08087,0],
        [0,0.08087]
        ])
        B=np.array([
        ])
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
        h=0.01
        mode='.'

class Verify:

    def verifyTemplate(At,Bt,h,mode,f,v,T,initset,uV,uP,strName):
        A=Benchmarks.createMatrix(At,Bt,h,mode)
        if (np.size(Bt)>0):
                n=np.size(At,1)+np.size(Bt,1)
        else:
                n=np.size(At,1)+0
        start_time = time.time()
        inp=InputFormation(A,n,f)
        '''print(A)
        print(f)
        exit(0)'''
        dyn=Dynamics(A,n,f)
        rg=StructLME([dyn.findConstStr()]+dyn.findVarBlkStr())
        if (rg.isCondSAT()):
            flg=True
            for i in range(1,int(T/h)+1):
                print("\n\n\n\n==============",i)
                dynRepA=DynamicsRepOp(dyn.findConstMat(),dyn.findVarBlkMats(),dyn.dimension)
                reach=CalcReachableSet(dynRepA,i,initset,f,v,uV,uP)
                flag=reach.isFeasible()
                end_time=time.time() - start_time
                if flag==True:
                    break
        else:
            flg=False

        print("\n\n------------ Final Report for ",strName," ------------")
        print("Satisfies Inf-Linear Conditions: ",flg)
        if (flg==True):
            print("Safety Status of the system upto time ",T,": ", end=" ")
            if (flag):
                print("Unsafe (Counter Example Given above)")
            else:
                print("Safe")
            print("Time Taken: ", end_time)
        else:
            print("Can't be verified!!")

    def verifyTemplateTimeVarying(At,Bt,h,mode,f,v,T,initset,uV,uP,strName,t):
        A=Benchmarks.createMatrix(At,Bt,h,mode)
        if (np.size(Bt)>0):
                n=np.size(At,1)+np.size(Bt,1)
        else:
                n=np.size(At,1)+0
        start_time = time.time()

        dyn=Dynamics(A,n,f)
        rg=StructLME([dyn.findConstStr()]+dyn.findVarBlkStr())
        if (rg.isCondSAT()):
            flg=True
            for i in range(1,int(T/h)+1):
                print("\n\n\n\n==============",i)
                dynRepA=DynamicsRepOp(dyn.findConstMat(),dyn.findVarBlkMats(),dyn.dimension)
                reach=CalcReachableSetTimeVarying(dynRepA,i,initset,f,v,uV,uP,t)
                flag=reach.isFeasible()
                end_time=time.time() - start_time
                if flag==True:
                    break
        else:
            flg=False

        print("\n\n------------ Final Report for ",strName," ------------")
        print("Satisfies Inf-Linear Conditions: ",flg)
        if (flg==True):
            print("Safety Status of the system upto time ",T,": ", end=" ")
            if (flag):
                print("Unsafe (Counter Example Given above)")
            else:
                print("Safe")
            print("Time Taken: ", end_time)
        else:
            print("Can't be verified!!")

    def verifyTemplateNonVerbose(At,Bt,h,mode,f,v,T,initset,uV,uP):

        report={
        "flg":False,
        "flag":False,
        "end_time":0
        }
        A=Benchmarks.createMatrix(At,Bt,h,mode)
        if (np.size(Bt)>0):
                n=np.size(At,1)+np.size(Bt,1)
        else:
                n=np.size(At,1)+0
        start_time = time.time()
        inp=InputFormation(A,n,f)
        dyn=Dynamics(A,n,f)
        rg=StructLME([dyn.findConstStr()]+dyn.findVarBlkStr())
        if (rg.isCondSAT()):
            flg=True
            for i in range(1,int(T/h)+1):
                dynRepA=DynamicsRepOp(dyn.findConstMat(),dyn.findVarBlkMats(),dyn.dimension)
                reach=CalcReachableSet(dynRepA,i,initset,f,v,uV,uP)
                flag=reach.isFeasible()
                end_time=time.time() - start_time
                if flag==True:
                    break
        else:
            flg=False

        report["flg"]=flg
        if flg:
            report["flag"]=flag
            report["end_time"]=end_time

        return report

    def verifyTemplateTimeVaryingNonVerbose(At,Bt,h,mode,f,v,T,initset,uV,uP,t):

        report={
        "flg":False,
        "flag":False,
        "end_time":0
        }
        A=Benchmarks.createMatrix(At,Bt,h,mode)
        if (np.size(Bt)>0):
                n=np.size(At,1)+np.size(Bt,1)
        else:
                n=np.size(At,1)+0
        start_time = time.time()
        inp=InputFormation(A,n,f)
        dyn=Dynamics(A,n,f)
        rg=StructLME([dyn.findConstStr()]+dyn.findVarBlkStr())
        if (rg.isCondSAT()):
            flg=True
            for i in range(1,int(T/h)+1):
                dynRepA=DynamicsRepOp(dyn.findConstMat(),dyn.findVarBlkMats(),dyn.dimension)
                reach=CalcReachableSetTimeVarying(dynRepA,i,initset,f,v,uV,uP,t)
                flag=reach.isFeasible()
                end_time=time.time() - start_time
                if flag==True:
                    break
        else:
            flg=False

        report["flg"]=flg
        if flg:
            report["flag"]=flag
            report["end_time"]=end_time

        return report

    def readFromFile(fname):
        '''A=
        B=
        h=
        mode=
        faults=
        variance=
        initialSet=
        strName=
        time=
        unsafeState=
        unsafeParam=
        Verify.verifyTemplate(A,B,h,mode,faults,variance,initialSet,strName,time,unsafeState,unsafeParam)'''
        tree = ET.parse(fname)
        root = tree.getroot()
        for child in root:
            #print(child.tag,child.text)
            if (child.tag=='h'):
                h=float(child.text.strip())
                print(h)
            elif (child.tag=='name'):
                strName=child.text.strip()
                print(strName)
            elif (child.tag=='mode'):
                mode=list(child.text.strip())[0]
                print(mode)
            elif (child.tag=='time'):
                time=int(child.text.strip())
                print(time)
            elif (child.tag=='unsafestate'):
                unsafeState=int(child.text.strip())
                print(unsafeState)
            elif (child.tag=='unsafeparam'):
                unsafeParam=float(child.text.strip())
                print(unsafeParam)
            elif (child.tag=='faults'):
                strFault=child.text.strip()
            elif (child.tag=='A'):
                strA=child.text.strip()
            elif (child.tag=='B'):
                strB=child.text.strip()
            elif (child.tag=='variance'):
                strVar=child.text.strip()
            elif (child.tag=='initialset'):
                strIS=child.text.strip()

        faults=[]
        for f in strFault.split('\n'):
            t=f.strip().split(',')
            faults.append((int(t[0].strip()),int(t[1].strip()),int(t[2].strip()),int(t[3].strip())))
        print(faults)

        variance=[]
        for v in strVar.split('\n'):
            t=v.strip().split(',')
            variance.append((int(t[0].strip()),int(t[1].strip())))
        print(variance)

        initialSet=[]
        for i in strIS.split('\n'):
            t=i.strip().split(',')
            initialSet.append((int(t[0].strip()),int(t[1].strip())))
        print(initialSet)

        a=[]
        for row in strA.split('\n'):
            tmp=row.strip().split(',')
            l=[]
            for t in tmp:
                if (len(t.strip())>0):
                    l.append(int(t.strip()))
            a.append(l)
        A=np.array(a)
        print(A)

        b=[]
        print(strB)
        for row in strB.split('\n'):
            tmp=row.strip().split(',')
            l=[]
            for t in tmp:
                if (len(t.strip())>0):
                    l.append(int(t.strip()))
            b.append(l)
        B=np.array(b)
        print(B)

        #print(int(time/h)+1)
        Verify.verifyTemplate(A,B,h,mode,faults,variance,time,initialSet,unsafeState,unsafeParam,strName)

    def readFromFileTimeVarying(fname):
        '''A=
        B=
        h=
        mode=
        faults=
        variance=
        initialSet=
        strName=
        time=
        unsafeState=
        unsafeParam=
        Verify.verifyTemplate(A,B,h,mode,faults,variance,initialSet,strName,time,unsafeState,unsafeParam)'''
        tree = ET.parse(fname)
        root = tree.getroot()
        for child in root:
            #print(child.tag,child.text)
            if (child.tag=='h'):
                h=float(child.text.strip())
                print(h)
            elif (child.tag=='name'):
                strName=child.text.strip()
                print(strName)
            elif (child.tag=='mode'):
                mode=list(child.text.strip())[0]
                print(mode)
            elif (child.tag=='time'):
                time=int(child.text.strip())
                print(time)
            elif (child.tag=='unsafestate'):
                unsafeState=int(child.text.strip())
                print(unsafeState)
            elif (child.tag=='unsafeparam'):
                unsafeParam=float(child.text.strip())
                print(unsafeParam)
            elif (child.tag=='faults'):
                strFault=child.text.strip()
            elif (child.tag=='A'):
                strA=child.text.strip()
            elif (child.tag=='B'):
                strB=child.text.strip()
            elif (child.tag=='variance'):
                strVar=child.text.strip()
            elif (child.tag=='initialset'):
                strIS=child.text.strip()
            elif (child.tag=='step'):
                step=int(child.text.strip())

        faults=[]
        for f in strFault.split('\n'):
            t=f.strip().split(',')
            faults.append((int(t[0].strip()),int(t[1].strip()),int(t[2].strip()),int(t[3].strip())))
        print(faults)

        variance=[]
        for v in strVar.split('\n'):
            l=[]
            for t in v.strip().split(';'):
                tmp=t.strip().split(',')
                l.append((int(tmp[0].strip()),int(tmp[1].strip())))
            variance.append(l)
        print(variance)

        initialSet=[]
        for i in strIS.split('\n'):
            t=i.strip().split(',')
            initialSet.append((int(t[0].strip()),int(t[1].strip())))
        print(initialSet)

        a=[]
        for row in strA.split('\n'):
            tmp=row.strip().split(',')
            l=[]
            for t in tmp:
                if (len(t.strip())>0):
                    l.append(int(t.strip()))
            a.append(l)
        A=np.array(a)
        print(A)

        b=[]
        print(strB)
        for row in strB.split('\n'):
            tmp=row.strip().split(',')
            l=[]
            for t in tmp:
                if (len(t.strip())>0):
                    l.append(int(t.strip()))
            b.append(l)
        B=np.array(b)
        print(B)

        #print(int(time/h)+1)
        Verify.verifyTemplateTimeVarying(A,B,h,mode,faults,variance,time,initialSet,unsafeState,unsafeParam,strName,step)

    def illustExample():
        faults=[(1,2,5,2),(3,2,5,2)]
        #faults=[(2,3,2,3)]
        variance=[(-0.01,0.01),(-0.02,0.02)]
        initialSet=[(1,2),(1,2),(3,4),(4,5),(4,5),(3,3),(2,2),(1,1)]
        strName="Illustrative Example Model"
        time=2
        unsafeState=1
        unsafeParam=77000
        Verify.verifyTemplate(Benchmarks.IllustExample.A,Benchmarks.IllustExample.B,Benchmarks.IllustExample.h,Benchmarks.IllustExample.mode,faults,variance,time,initialSet,unsafeState,unsafeParam,strName)

    def illustExample2():
        faults=[(0,1,1,1),(0,1,2,1)]
        #faults=[(2,3,2,3)]
        variance=[(-0.01,0.01),(-0.02,0.02)]
        initialSet=[(1,2),(2,2),(3,3),(1,1)]
        strName="Illustrative Example Model"
        time=2
        unsafeState=0
        unsafeParam=77000
        Verify.verifyTemplate(Benchmarks.IllustExample2.A,Benchmarks.IllustExample2.B,Benchmarks.IllustExample2.h,Benchmarks.IllustExample2.mode,faults,variance,time,initialSet,unsafeState,unsafeParam,strName)

    def test():
        faults=[(3,3,0,3)]
        #faults=[(2,3,2,3)]
        variance=[(5,5)]
        initialSet=[(-1,1),(-1,1),(-1,1),(-1,1),(-1,1),(-1,1)]
        strName="Example Model"
        time=20
        unsafeState=0
        unsafeParam=0
        Verify.verifyTemplate(Benchmarks.Test.A,Benchmarks.Test.B,Benchmarks.Test.h,Benchmarks.Test.mode,faults,variance,time,initialSet,unsafeState,unsafeParam,strName)

    def testTimevarying():
        faults=[(3,3,0,3)]
        #faults=[(2,3,2,3)]
        variance=[[(5,5),(1,1)]]
        initialSet=[(-1,1),(-1,1),(-1,1),(-1,1),(-1,1),(-1,1)]
        strName="Example Model (Time Varying)"
        time=20
        step=1000
        unsafeState=0
        unsafeParam=0
        Verify.verifyTemplateTimeVarying(Benchmarks.Test.A,Benchmarks.Test.B,Benchmarks.Test.h,Benchmarks.Test.mode,faults,variance,time,initialSet,unsafeState,unsafeParam,strName,step)

    def PKPD2():
        faults=[(0,4,4,1)]
        #faults=[(3,1,4,1)]
        variance=[(5,5)]
        #initialSet=[(1,6),(0,10),(0,10),(1,8),(0,200)]
        initialSet=[(1,6),(0,10),(0,10),(1,8),(1,1)]
        strName="PK/PD Model - II"
        time=20
        unsafeState=0
        unsafeParam=0
        Verify.verifyTemplate(Benchmarks.PKPD2.A,Benchmarks.PKPD2.B,Benchmarks.PKPD2.h,Benchmarks.PKPD2.mode,faults,variance,time,initialSet,unsafeState,unsafeParam,strName)

    def PKPD2TimeVarying():
        faults=[(0,4,4,1)]
        #faults=[(3,1,4,1)]
        #variance=[[(-1000,1002),(1,1),(0.1,1.9),(-1000,1002),(-1000,1002),(-1,2),(-2,3),(-100,102),(0.8,1.2),(0.5,1.5),(-1,2),(0.1,1.9),(-1,2),(0.5,1.5),(0.6,0.4),(0.8,1.2),(0.7,1.3),(-1,2),(0.8,1.2),(0.5,1.5)]]
        variance=[[(-1000,2000),(-1000,2000),(-1000,2000),(-1000,2000),(-1000,2000),(-1000,2000),(-1000,2000),(-1000,2000),(-1000,2000),(-1000,2000),(-1000,2000),(-1000,2000),(-1000,2000),(-1000,2000),(-1000,2000),(-1000,2000),(-1000,2000),(-1000,2000),(-1000,2000),(-1000,2000),(-1000,2000)]]
        initialSet=[(1,6),(0,10),(0,10),(1,8),(1,1)]
        var=30
        strName="PK/PD Model - II (Time-Varying Perturbation)"
        time=20
        unsafeState=0
        unsafeParam=0
        Verify.verifyTemplateTimeVarying(Benchmarks.PKPD2.A,Benchmarks.PKPD2.B,Benchmarks.PKPD2.h,Benchmarks.PKPD2.mode,faults,variance,time,initialSet,unsafeState,unsafeParam,strName,var)

    def flightEnvelope():
        '''faults=[(0,7,10,6)]
        #faults=[(0,8,8,8)]
        variance=[(5,5)]
        #initialSet=[(3,4),(3,4),(3,5),(3,4),(3,4),(3,5),(-0.50,0.50),(-0.50,0.50),(-0.50,0.50),(-1,1),(-1,1),(-1,1),(5,10),(-0.1,0.1),(-0.1,0.1),(-0.1,0.1)]
        #initialSet=[(3,4),(3,4),(3,5),(3,4),(3,4),(3,5),(-0.50,0.50),(-0.50,0.50),(-0.50,0.50),(-1,1),(-1,1),(-1,1),(1,1),(1,1),(1,1),(1,1)]
        initialSet=[(3,4),(3,4),(3,5),(3,4),(3,4),(3,5),(-0.50,0.50),(-0.50,0.50),(-0.50,0.50),(-1,1),(1,1),(1,1),(1,1),(1,1),(1,1),(1,1)]
        strName="Flight Envelope Model"
        time=20
        unsafeState=0
        unsafeParam=0
        Verify.verifyTemplate(Benchmarks.FlightEnvelope.A,Benchmarks.FlightEnvelope.B,Benchmarks.FlightEnvelope.h,Benchmarks.FlightEnvelope.mode,faults,variance,time,initialSet,unsafeState,unsafeParam,strName)
        '''
        faults=[(0,6,6,10)]
        #faults=[(0,8,8,8)]
        #faults=[(0,2,1,14)]
        variance=[(1,1)]
        #initialSet=[(3,4),(3,4),(3,5),(3,4),(3,4),(3,5),(-0.50,0.50),(-0.50,0.50),(-0.50,0.50),(-1,1),(-1,1),(-1,1),(5,10),(-0.1,0.1),(-0.1,0.1),(-0.1,0.1)]
        #initialSet=[(3,4),(3,4),(3,5),(3,4),(3,4),(3,5),(-0.50,0.50),(-0.50,0.50),(-0.50,0.50),(-1,1),(-1,1),(-1,1),(1,1),(1,1),(1,1),(1,1)]
        #initialSet=[(3,4),(3,4),(3,5),(3,4),(3,4),(3,5),(-0.50,0.50),(-0.50,0.50),(-0.50,0.50),(-1,1),(1,1),(1,1),(1,1),(1,1),(1,1),(1,1)]
        #initialSet=[(3,4),(3,4),(3,4),(3,4),(3,4),(3,5),(-0.50,0.50),(1,1),(1,1),(1,1),(1,1),(1,1),(1,1),(1,1),(1,1),(1,1)]
        initialSet=[(3,4),(3,4),(3,4),(3,4),(3,4),(3,5),(-0.5,0.5),(1,1),(1,1),(1,1),(1,1),(1,1),(1,1),(1,1),(1,1),(1,1)]
        strName="Flight Envelope Model"
        time=20
        unsafeState=3
        unsafeParam=-15
        Verify.verifyTemplate(Benchmarks.FlightEnvelope.A,Benchmarks.FlightEnvelope.B,Benchmarks.FlightEnvelope.h,Benchmarks.FlightEnvelope.mode,faults,variance,time,initialSet,unsafeState,unsafeParam,strName)

    def flightEnvelopeTimeVarying():
        faults=[(0,6,6,10)]
        #faults=[(0,8,8,8)]
        #faults=[(0,2,1,14)]
        #variance=[[(0.99,1.01),(0.98,1.02),(1,1),(0.9,1.1),(0.95,1.05),(0.8,1.2),(0.7,1.3),(0.9,1.1),(1,1),(0.99,1.01),(0.95,1.05),(1,1),(0.98,1.02),(0.8,1.2),(1,1),(0.75,1.25),(1,1),(0.9,1.1),(0.7,1.3),(0.9,1.1)]]
        variance=[[(0.87,1.13),(0.87,1.13),(0.87,1.13),(0.87,1.13),(0.87,1.13),(0.87,1.13),(0.87,1.13),(0.87,1.13),(0.87,1.13),(0.87,1.13),(0.87,1.13),(0.87,1.13),(0.87,1.13),(0.87,1.13),(0.87,1.13),(0.87,1.13),(0.87,1.13),(0.87,1.13),(0.87,1.13),(0.87,1.13)]]
        initialSet=[(3,4),(3,4),(3,4),(3,4),(3,4),(3,5),(-0.50,0.50),(1,1),(1,1),(1,1),(1,1),(1,1),(1,1),(1,1),(1,1),(1,1)]
        strName="Flight Envelope Model (Time-Varying Perturbation)"
        time=20
        var=6
        unsafeState=3
        unsafeParam=-15

        Verify.verifyTemplateTimeVarying(Benchmarks.FlightEnvelope.A,Benchmarks.FlightEnvelope.B,Benchmarks.FlightEnvelope.h,Benchmarks.FlightEnvelope.mode,faults,variance,time,initialSet,unsafeState,unsafeParam,strName,var)

    def coOPVehiclesI():
        faults=[(0,9,9,1)]
        variance=[(5,5)]
        #initialSet=[(0,0),(0,0),(0,0),(0,0),(0,0),(0,0),(0,0),(0,0),(0,0),(1,1)]
        initialSet=[(1,1),(1,1),(1,1),(1,1),(1,1),(1,1),(1,1),(1,1),(1,1),(1,1)]
        strName="Networked Cooperative Platoon of Vehicles - I"
        time=20
        unsafeState=0
        unsafeParam=-26
        Verify.verifyTemplate(Benchmarks.CoOPVehiclesI.A,Benchmarks.CoOPVehiclesI.B,Benchmarks.CoOPVehiclesI.h,Benchmarks.CoOPVehiclesI.mode,faults,variance,time,initialSet,unsafeState,unsafeParam,strName)

    def coOPVehiclesII():
        faults=[(0,9,9,1)]
        variance=[(5,5)]
        initialSet=[(1,1),(1,1),(1,1),(1,1),(1,1),(1,1),(1,1),(1,1),(1,1),(1,1)]
        strName="Networked Cooperative Platoon of Vehicles - II"
        time=20
        unsafeState=0
        unsafeParam=0
        Verify.verifyTemplate(Benchmarks.CoOPVehiclesII.A,Benchmarks.CoOPVehiclesII.B,Benchmarks.CoOPVehiclesII.h,Benchmarks.CoOPVehiclesII.mode,faults,variance,time,initialSet,unsafeState,unsafeParam,strName)

    def coOPVehiclesIITimeVarying():
        faults=[(0,9,9,1)]
        #variance=[[(0.5,1.5),(0.6,1.4),(0.5,1.5),(0.7,1.3),(-1,2),(0.5,1.5),(0.4,0.6),(0.5,1.5),(0.6,1.4),(0.7,1.3),(0.65,1.35),(0.6,1.4),(0.7,1.3),(0.6,1.4),(0.5,1.5),(0.7,1.3),(-1,2),(0.7,1.3),(0.5,1.5),(0.8,1.2)]]
        variance=[[(0.1,1.9),(0.1,1.9),(0.1,1.9),(0.1,1.9),(0.1,1.9),(0.1,1.9),(0.1,1.9),(0.1,1.9),(0.1,1.9),(0.1,1.9),(0.1,1.9),(0.1,1.9),(0.1,1.9),(0.1,1.9),(0.1,1.9),(0.1,1.9),(0.1,1.9),(0.1,1.9),(0.1,1.9),(0.1,1.9)]]
        initialSet=[(1,1),(1,1),(1,1),(1,1),(1,1),(1,1),(1,1),(1,1),(1,1),(1,1)]
        strName="Networked Cooperative Platoon of Vehicles - II (Time-Varying Perturbation)"
        time=20
        var=20
        unsafeState=0
        unsafeParam=1
        Verify.verifyTemplateTimeVarying(Benchmarks.CoOPVehiclesII.A,Benchmarks.CoOPVehiclesII.B,Benchmarks.CoOPVehiclesII.h,Benchmarks.CoOPVehiclesII.mode,faults,variance,time,initialSet,unsafeState,unsafeParam,strName,var)

    def holes():
        faults=[(0, 4, 4, 4)]
        variance = [(5,5)]
        #initialSet=[(0,0),(0,0),(0,0),(0,0),(0,0),(-100,100),(0,0),(0,0),(-100,100),(-100,100)]
        initialSet=[(0,0),(0,0),(0,0),(0,0),(0,0),(0,0),(-100,-100),(-100,-100)]
        strName="Holes Plant Xp"
        time=20
        unsafeState=0
        unsafeParam=0
        Verify.verifyTemplate(Benchmarks.HolesPXp.A,Benchmarks.HolesPXp.B,Benchmarks.HolesPXp.h,Benchmarks.HolesPXp.mode,faults,variance,time,initialSet,unsafeState,unsafeParam,strName)

    def motorTransmission1():
        faults=[(0,4,4,3)]
        variance = [(5,5)]
        #initialSet=[(0,0),(-0.08,0.08),(-0.165,-0.165),(-0.01,0.01),(0,0),(70,70),(1,1)]
        initialSet=[(0,0),(-0.08,0.08),(-0.0165,-0.0165),(-0.01,0.01),(0,0),(70,70),(1,1)]
        strName="Motor Transmission - I"
        time=20
        unsafeState=0
        unsafeParam=0
        Verify.verifyTemplate(Benchmarks.MotorTransmission1.A,Benchmarks.MotorTransmission1.B,Benchmarks.MotorTransmission1.h,Benchmarks.MotorTransmission1.mode,faults,variance,time,initialSet,unsafeState,unsafeParam,strName)

    def motorTransmission1TimeVarying():
        faults=[(0,4,4,3)]
        #variance = [[(0.7,1.3),(-0.2,1),(0.01,1.99),(-0.1,1),(0.8,1.2),(-0.3,1),(0.7,1.3),(1,1),(0.8,1.2),(0.6,0.4),(0.65,0.35),(0.8,1.2),(1,1),(1,1.2),(0.7,1),(0.8,1),(1,1.1),(1,1.2),(1,1),(0.9,1.1)]]
        variance = [[(-0.01,1),(-0.01,1),(-0.01,1),(-0.01,1),(-0.01,1),(-0.01,1),(-0.01,1),(-0.01,1),(-0.01,1),(-0.01,1),(-0.01,1),(-0.01,1),(-0.01,1),(-0.01,1),(-0.01,1),(-0.01,1),(-0.01,1),(-0.01,1),(-0.01,1),(-0.01,1)]]
        #initialSet=[(0,0),(-0.08,0.08),(-0.165,-0.165),(-0.01,0.01),(0,0),(70,70),(1,1)]
        initialSet=[(0,0),(-0.08,0.08),(-0.0165,-0.0165),(-0.01,0.01),(0,0),(70,70),(1,1)]
        var=2
        strName="Motor Transmission - I (Time-Varying Perturbation)"
        time=20
        unsafeState=2
        unsafeParam=-0.02
        Verify.verifyTemplateTimeVarying(Benchmarks.MotorTransmission1.A,Benchmarks.MotorTransmission1.B,Benchmarks.MotorTransmission1.h,Benchmarks.MotorTransmission1.mode,faults,variance,time,initialSet,unsafeState,unsafeParam,strName,var)

    def motorTransmission2():
        faults=[(2,3,0,2)]
        variance = [(5,5)]
        initialSet=[(0,0),(-0.08,0.08),(-0.165,-0.165),(-0.01,0.01),(0,0)]
        strName="Motor Transmission - II"
        time=20
        unsafeState=0
        unsafeParam=0
        Verify.verifyTemplate(Benchmarks.MotorTransmission2.A,Benchmarks.MotorTransmission2.B,Benchmarks.MotorTransmission2.h,Benchmarks.MotorTransmission2.mode,faults,variance,time,initialSet,unsafeState,unsafeParam,strName)

class VerifyFaultFree:

    def verifyTemplateFaultFree(At,Bt,h,mode,T,initset,unsafeV,unsafeP,strName):
        A=Benchmarks.createMatrix(At,Bt,h,mode)
        if (np.size(Bt)>0):
                n=np.size(At,1)+np.size(Bt,1)
        else:
                n=np.size(At,1)+0
        start_time = time.time()
        for i in range(1,int(T/h)+1):
            print("\n\n\n\n==============",i)
            reach=CalcReachableSetFaultFree(A,n,i,initset,unsafeV,unsafeP)
            flag=reach.isFeasible()
            end_time=time.time() - start_time
            if flag==True:
                break
            else:
                flg=False

        print("\n\n------------ Final Report for ",strName,"(Without Any Fault) ------------")
        print("Safety Status of the system upto time ",T,": ", end=" ")
        if (flag):
            print("Unsafe (Counter Example Given above)")
        else:
            print("Safe")
        print("Time Taken: ", end_time)

    def verifyTemplateFaultFreeNonVerbose(At,Bt,h,mode,T,initset,unsafeV,unsafeP):
        report={
        "flag":False,
        "end_time":0
        }
        A=Benchmarks.createMatrix(At,Bt,h,mode)
        if (np.size(Bt)>0):
                n=np.size(At,1)+np.size(Bt,1)
        else:
                n=np.size(At,1)+0
        start_time = time.time()
        flg=True
        for i in range(1,int(T/h)+1):
            print("\n\n\n\n==============",i)
            reach=CalcReachableSetFaultFree(A,n,i,initset,unsafeV,unsafeP)
            flag=reach.isFeasible()
            end_time=time.time() - start_time
            if flag==True:
                break
            else:
                flg=False

        report["flag"]=flag
        report["end_time"]=end_time

        return report

    def readFromFile(fname):
        tree = ET.parse(fname)
        root = tree.getroot()
        for child in root:
            #print(child.tag,child.text)
            if (child.tag=='h'):
                h=float(child.text.strip())
                print(h)
            elif (child.tag=='name'):
                strName=child.text.strip()
                print(strName)
            elif (child.tag=='mode'):
                mode=list(child.text.strip())[0]
                print(mode)
            elif (child.tag=='time'):
                time=int(child.text.strip())
                print(time)
            elif (child.tag=='unsafestate'):
                unsafeState=int(child.text.strip())
                print(unsafeState)
            elif (child.tag=='unsafeparam'):
                unsafeParam=float(child.text.strip())
                print(unsafeParam)
            elif (child.tag=='A'):
                strA=child.text.strip()
            elif (child.tag=='B'):
                strB=child.text.strip()
            elif (child.tag=='initialset'):
                strIS=child.text.strip()

        initialSet=[]
        for i in strIS.split('\n'):
            t=i.strip().split(',')
            initialSet.append((int(t[0].strip()),int(t[1].strip())))
        print(initialSet)

        a=[]
        for row in strA.split('\n'):
            tmp=row.strip().split(',')
            l=[]
            for t in tmp:
                if (len(t.strip())>0):
                    l.append(int(t.strip()))
            a.append(l)
        A=np.array(a)
        print(A)

        b=[]
        print(strB)
        for row in strB.split('\n'):
            tmp=row.strip().split(',')
            l=[]
            for t in tmp:
                if (len(t.strip())>0):
                    l.append(int(t.strip()))
            b.append(l)
        B=np.array(b)
        print(B)

        VerifyFaultFree.verifyTemplateFaultFree(A,B,h,mode,time,initialSet,unsafeState,unsafeParam,strName)

    def illustExample():
        initialSet=[(-1,1),(-1,1),(-1,1),(-1,1),(-1,1),(-1,1)]
        strName="Illustrative Example Model"
        time=5
        uV=0
        uP=0
        VerifyFaultFree.verifyTemplateFaultFree(Benchmarks.IllustExample.A,Benchmarks.IllustExample.B,Benchmarks.IllustExample.h,Benchmarks.IllustExample.mode,time,initialSet,uV,uP,strName)

    def test():
        initialSet=[(-1,1),(-1,1),(-1,1),(-1,1),(-1,1),(-1,1)]
        strName="Example Model"
        time=20
        uV=0
        uP=0
        VerifyFaultFree.verifyTemplateFaultFree(Benchmarks.Test.A,Benchmarks.Test.B,Benchmarks.Test.h,Benchmarks.Test.mode,time,initialSet,uV,uP,strName)

    def PKPD2():
        initialSet=[(1,6),(0,10),(0,10),(1,8),(0,200)]
        strName="PK/PD Model - II"
        time=20
        uV=0
        uP=0
        VerifyFaultFree.verifyTemplateFaultFree(Benchmarks.PKPD2.A,Benchmarks.PKPD2.B,Benchmarks.PKPD2.h,Benchmarks.PKPD2.mode,time,initialSet,uV,uP,strName)
        #VerifyFaultFree.verifyTemplate(Benchmarks.PKPD2.A,Benchmarks.PKPD2.B,Benchmarks.PKPD2.n1,Benchmarks.PKPD2.n2,Benchmarks.PKPD2.h,Benchmarks.PKPD2.mode,faults,variance,time,initialSet,[],[],0,strName)

    def flightEnvelope():
        initialSet=[(3,4),(3,4),(3,5),(3,4),(3,4),(3,5),(-0.50,0.50),(-0.50,0.50),(-0.50,0.50),(-1,1),(1,1),(1,1),(1,1),(1,1),(1,1),(1,1)]
        strName="Flight Envelope Model"
        time=20
        uV=0
        uP=0
        VerifyFaultFree.verifyTemplateFaultFree(Benchmarks.FlightEnvelope.A,Benchmarks.FlightEnvelope.B,Benchmarks.FlightEnvelope.h,Benchmarks.FlightEnvelope.mode,time,initialSet,uV,uP,strName)
        #VerifyFaultFree.verifyTemplate(Benchmarks.FlightEnvelope.A,Benchmarks.FlightEnvelope.B,Benchmarks.FlightEnvelope.n1,Benchmarks.FlightEnvelope.n2,Benchmarks.FlightEnvelope.h,Benchmarks.FlightEnvelope.mode,faults,variance,time,initialSet,[],[],0,strName)

    def coOPVehiclesI():
        initialSet=[(1,1),(1,1),(1,1),(1,1),(1,1),(1,1),(1,1),(1,1),(1,1),(-9,1)]
        strName="Networked Cooperative Platoon of Vehicles - I"
        time=20
        uV=0
        uP=0
        VerifyFaultFree.verifyTemplateFaultFree(Benchmarks.CoOPVehiclesI.A,Benchmarks.CoOPVehiclesI.B,Benchmarks.CoOPVehiclesI.h,Benchmarks.CoOPVehiclesI.mode,time,initialSet,uV,uP,strName)
        #Benchmarks.verifyTemplateFaultFree(Benchmarks.CoOPVehiclesI.A,Benchmarks.CoOPVehiclesI.B,Benchmarks.CoOPVehiclesI.n1,Benchmarks.CoOPVehiclesI.n2,Benchmarks.CoOPVehiclesI.h,Benchmarks.CoOPVehiclesI.mode,time,initialSet,strName)
        #VerifyFaultFree.verifyTemplate(Benchmarks.CoOPVehiclesI.A,Benchmarks.CoOPVehiclesI.B,Benchmarks.CoOPVehiclesI.n1,Benchmarks.CoOPVehiclesI.n2,Benchmarks.CoOPVehiclesI.h,Benchmarks.CoOPVehiclesI.mode,faults,variance,time,initialSet,[],[],0,strName)

    def coOPVehiclesII():
        strName="Networked Cooperative Platoon of Vehicles - II"
        time=20
        uV=0
        uP=0
        VerifyFaultFree.verifyTemplateFaultFree(Benchmarks.CoOPVehiclesII.A,Benchmarks.CoOPVehiclesII.B,Benchmarks.CoOPVehiclesII.h,Benchmarks.CoOPVehiclesII.mode,time,initialSet,uV,uP,strName)
        #VerifyFaultFree.verifyTemplate(Benchmarks.CoOPVehiclesII.A,Benchmarks.CoOPVehiclesII.B,Benchmarks.CoOPVehiclesII.n1,Benchmarks.CoOBenchmarks.CoOPVehiclesI.n1,Benchmarks.CoOPVehiclesI.n2PVehiclesII.n2,Benchmarks.CoOPVehiclesII.h,Benchmarks.CoOPVehiclesII.mode,faults,variance,time,initialSet,[],[],0,strName)

    def holes():
        initialSet=[(0,0),(0,0),(0,0),(0,0),(0,0),(0,0),(-100,100),(-100,100)]
        strName="Holes Plant Xp"
        time=20
        uV=0
        uP=0
        VerifyFaultFree.verifyTemplateFaultFree(Benchmarks.HolesPXp.A,Benchmarks.HolesPXp.B,Benchmarks.HolesPXp.h,Benchmarks.HolesPXp.mode,time,initialSet,uV,uP,strName)
        #VerifyFaultFree.verifyTemplate(Benchmarks.HolesPXp.A,Benchmarks.HolesPXp.B,Benchmarks.HolesPXp.n1,Benchmarks.HolesPXp.n2,Benchmarks.HolesPXp.h,Benchmarks.HolesPXp.mode,faults,variance,time,initialSet,[],[],0,strName)

    def motorTransmission1():
        initialSet=[(0,0),(-0.08,0.08),(-0.0165,-0.0165),(-0.01,0.01),(0,0),(70,70),(1,1)]
        strName="Motor Transmission - I"
        time=20
        uV=0
        uP=0
        VerifyFaultFree.verifyTemplateFaultFree(Benchmarks.MotorTransmission1.A,Benchmarks.MotorTransmission1.B,Benchmarks.MotorTransmission1.h,Benchmarks.MotorTransmission1.mode,time,initialSet,uV,uP,strName)
        #VerifyFaultFree.verifyTemplate(Benchmarks.MotorTransmission1.A,Benchmarks.MotorTransmission1.B,Benchmarks.MotorTransmission1.n1,Benchmarks.MotorTransmission1.n2,Benchmarks.MotorTransmission1.h,Benchmarks.MotorTransmission1.mode,faults,variance,time,initialSet,[],[],0,strName)

    def motorTransmission2():
        initialSet=[(0,0),(-0.08,0.08),(-0.165,-0.165),(-0.01,0.01),(0,0)]
        strName="Motor Transmission - II"
        time=20
        uV=0
        uP=0
        VerifyFaultFree.verifyTemplateFaultFree(Benchmarks.MotorTransmission2.A,Benchmarks.MotorTransmission2.B,Benchmarks.MotorTransmission2.h,Benchmarks.MotorTransmission2.mode,faults,variance,time,initialSet,uV,uP,strName)

class Finder:
    def __init__(self,A,B,h,mode):
        self.A=A
        self.B=B
        self.appxM=self.createMatrix(A,B,h,mode)
        if (np.size(B)>1):
            self.n=np.size(A,1)+np.size(B,1)
        else:
            self.n=np.size(A,1)+0

    def createMatrix(self,A,B,h,mode):
        ''' Creates a single matrix based on
        . or +.
        In case of . a roungh approximation is
        done'''

        n1=np.size(A,1)
        if (np.size(B)>0):
            n2=np.size(B,1)
        else:
            n2=0
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
                    print("\n\n----",m,k,"(Block Size) ----")
                    print(flt)
                    print("====================\n")
                else:
                    print("Finished ",m,",",k)

    def viewStructLME(self,flt):
        inp=InputFormation(self.appxM,self.n,flt)
        dyn=Dynamics(inp.formMatA(),self.n,flt)
        rg=StructLME([dyn.findConstStr()]+dyn.findVarBlkStr())
        print(rg.isCondSAT())
        rg.visualizeStructLME()

    def enumBlock(self,m,k):
        '''Enumerates possible blocks of
        size m*k'''
        l=m
        b=k
        fltList=[]
        for i in range(self.n-l+1):
            for j in range(self.n-b+1):
                flt=(i,l,j,b)
                inp=InputFormation(self.appxM,self.n,[flt])
                dyn=Dynamics(inp.formMatA(),self.n,[flt])
                #print(dyn.findConstStr())
                rg=StructLME([dyn.findConstStr()]+dyn.findVarBlkStr())
                if (rg.isCondSAT()):
                    fltList.append(flt)
                    '''print(flt)
                    rg.visualizeStructLME()
                    print('\n')'''
        fltList=Finder.filterOutInt(fltList)
        #fltList=[(1,2,3,2),(4,5,6,7),(8,9,10,11),(3,4,5,6)]

        if len(fltList)>0:
            for n in range(len(fltList),0,-1):
                for elem in itertools.combinations(fltList,n):
                    #print(elem)
                    inp=InputFormation(self.appxM,self.n,elem)
                    dyn=Dynamics(inp.formMatA(),self.n,elem)
                    rg=StructLME([dyn.findConstStr()]+dyn.findVarBlkStr())
                    if (rg.isCondSAT()):
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

class Search:

    def searchTemplate(A,B,h,mode,flt=None):
        fnd=Finder(A,B,h,mode)
        print("-----------")
        if flt==None:
            start_time=time.time()
            fnd.enumAllBlocks()
            print("==========")
            print("Time Taken: ",time.time()-start_time)
        else:
            fnd.viewStructLME(flt)

    def readFromFile(fname):
        tree = ET.parse(fname)
        root = tree.getroot()
        for child in root:
            if (child.tag=='h'):
                h=float(child.text.strip())
                print(h)
            elif (child.tag=='mode'):
                mode=list(child.text.strip())[0]
                print(mode)
            elif (child.tag=='faults'):
                strFault=child.text.strip()
            elif (child.tag=='A'):
                strA=child.text.strip()
            elif (child.tag=='B'):
                strB=child.text.strip()

        faults=[]
        for f in strFault.split('\n'):
            t=f.strip().split(',')
            #print(t)
            if (len(t[0].strip())>0):
                    faults.append((int(t[0].strip()),int(t[1].strip()),int(t[2].strip()),int(t[3].strip())))
        print(faults)

        a=[]
        for row in strA.split('\n'):
            tmp=row.strip().split(',')
            l=[]
            for t in tmp:
                if (len(t.strip())>0):
                    l.append(int(t.strip()))
            a.append(l)
        A=np.array(a)
        print(A)

        b=[]
        #print(strB)
        for row in strB.split('\n'):
            tmp=row.strip().split(',')
            l=[]
            for t in tmp:
                if (len(t.strip())>0):
                    l.append(int(t.strip()))
            b.append(l)
        B=np.array(b)
        print(B)

        if (faults==[]):
            faults=None
        Search.searchTemplate(A,B,h,mode,faults)

    def illustExample(flt=None):
        Search.searchTemplate(Benchmarks.IllustExample.A,Benchmarks.IllustExample.B,Benchmarks.IllustExample.h,Benchmarks.IllustExample.mode,flt)

    def illustExample2(flt=None):
        Search.searchTemplate(Benchmarks.IllustExample2.A,Benchmarks.IllustExample2.B,Benchmarks.IllustExample2.h,Benchmarks.IllustExample2.mode,flt)

    def madeUp(flt=None):
        Search.searchTemplate(Benchmarks.Test.A,Benchmarks.Test.B,Benchmarks.Test.h,Benchmarks.Test.mode,flt)

    def flightEnvelope(flt=None):
        Search.searchTemplate(Benchmarks.FlightEnvelope.A,Benchmarks.FlightEnvelope.B,Benchmarks.FlightEnvelope.h,Benchmarks.FlightEnvelope.mode,flt)

    def dcConverter(flt=None):
        Search.searchTemplate(Benchmarks.DCConv.A,Benchmarks.DCConv.B,Benchmarks.DCConv.h,Benchmarks.DCConv.mode,flt)

    def fiveVehiclePlatoon(flt=None):
        Search.searchTemplate(Benchmarks.FiveVehiclePlatton.A,Benchmarks.FiveVehiclePlatton.B,Benchmarks.FiveVehiclePlatton.h,Benchmarks.FiveVehiclePlatton.mode,flt)

    def tenVehiclePlatoon(flt=None):
        Search.searchTemplate(Benchmarks.TenVehiclePlatton.A,Benchmarks.TenVehiclePlatton.B,Benchmarks.TenVehiclePlatton.h,Benchmarks.TenVehiclePlatton.mode,flt)

    def coOPVehicles(flt=None):
        Search.searchTemplate(Benchmarks.CoOPVehiclesI.A,Benchmarks.CoOPVehiclesI.B,Benchmarks.CoOPVehiclesI.h,Benchmarks.CoOPVehiclesI.mode,flt)

    def pkpd(flt=None):
        Search.searchTemplate(Benchmarks.PKPD.A,Benchmarks.PKPD.B,Benchmarks.PKPD.h,Benchmarks.PKPD.mode,flt)

    def pkpd2(flt=None):
        Search.searchTemplate(Benchmarks.PKPD2.A,Benchmarks.PKPD2.B,Benchmarks.PKPD2.h,Benchmarks.PKPD2.mode,flt)

    def spaceCraftRndzvs(flt=None):
        Search.searchTemplate(Benchmarks.SpaceCraftRndzvs.A,Benchmarks.SpaceCraftRndzvs.B,Benchmarks.SpaceCraftRndzvs.h,Benchmarks.SpaceCraftRndzvs.mode,flt)

    def holesCXc(flt=None):
        Search.searchTemplate(Benchmarks.HolesCXc.A,Benchmarks.HolesCXc.B,Benchmarks.HolesCXc.h,Benchmarks.HolesCXc.mode,flt)

    def holesPDp(flt=None):
        Search.searchTemplate(Benchmarks.HolesPDp.A,Benchmarks.HolesPDp.B,Benchmarks.HolesPDp.h,Benchmarks.HolesPDp.mode,flt)

    def holesPXp(flt=None):
        Search.searchTemplate(Benchmarks.HolesPXp.A,Benchmarks.HolesPXp.B,Benchmarks.HolesPXp.h,Benchmarks.HolesPXp.mode,flt)

    def motorTransmission1(flt=None):
        Search.searchTemplate(Benchmarks.MotorTransmission1.A,Benchmarks.MotorTransmission1.B,Benchmarks.MotorTransmission1.h,Benchmarks.MotorTransmission1.mode,flt)

    def motorTransmission2(flt=None):
        Search.searchTemplate(Benchmarks.MotorTransmission2.A,Benchmarks.MotorTransmission2.B,Benchmarks.MotorTransmission2.h,Benchmarks.MotorTransmission2.mode,flt)

class ConsolidatedVerificationReport:

    def conVerTemplate(At,Bt,h,mode,f,v,T,initset,uV,uP,strName):
        reportWithFault=Verify.verifyTemplateNonVerbose(At,Bt,h,mode,f,v,T,initset,uV,uP)
        reportWithouFault=VerifyFaultFree.verifyTemplateFaultFreeNonVerbose(At,Bt,h,mode,T,initset,uV,uP)

        print("\n\n------------ Final Report for ",strName,"(With Faults) ------------")
        print("Satisfies Inf-Linear Conditions: ",reportWithFault["flg"])
        if (reportWithFault["flg"]==True):
            print("Safety Status of the system upto time ",T,": ", end=" ")
            if (reportWithFault["flag"]):
                print("Unsafe (To get the Counter Example please use Verify.verifyTemplate(...))")
            else:
                print("Safe")
            print("Time Taken: ", reportWithFault["end_time"])
        else:
            print("Can't be verified!!")

        print("\n\n------------ Final Report for ",strName,"(Without Any Fault) ------------")
        print("Safety Status of the system upto time ",T,": ", end=" ")
        if (reportWithouFault["flag"]):
            print("Unsafe (To get the Counter Example use VerifyFaultFree.verifyTemplate(...))")
        else:
            print("Safe")
        print("Time Taken: ", reportWithouFault["end_time"])

    def conVerTemplateTimeVarying(At,Bt,h,mode,f,v,T,initset,uV,uP,strName,t):
        reportWithFault=Verify.verifyTemplateTimeVaryingNonVerbose(At,Bt,h,mode,f,v,T,initset,uV,uP,t)
        reportWithouFault=VerifyFaultFree.verifyTemplateFaultFreeNonVerbose(At,Bt,h,mode,T,initset,uV,uP)

        print("\n\n------------ Final Report for ",strName,"(With Time-Varying Faults) ------------")
        print("Satisfies Inf-Linear Conditions: ",reportWithFault["flg"])
        if (reportWithFault["flg"]==True):
            print("Safety Status of the system upto time ",T,": ", end=" ")
            if (reportWithFault["flag"]):
                print("Unsafe (To get the Counter Example please use Verify.verifyTemplateTimeVarying(...))")
            else:
                print("Safe")
            print("Time Taken: ", reportWithFault["end_time"])
        else:
            print("Can't be verified!!")

        print("\n\n------------ Final Report for ",strName,"(Without Any Fault) ------------")
        print("Safety Status of the system upto time ",T,": ", end=" ")
        if (reportWithouFault["flag"]):
            print("Unsafe (To get the Counter Example use VerifyFaultFree.verifyTemplate(...))")
        else:
            print("Safe")
        print("Time Taken: ", reportWithouFault["end_time"])

    def readFromFile(fname):
        tree = ET.parse(fname)
        root = tree.getroot()
        for child in root:
            #print(child.tag,child.text)
            if (child.tag=='h'):
                h=float(child.text.strip())
                print(h)
            elif (child.tag=='name'):
                strName=child.text.strip()
                print(strName)
            elif (child.tag=='mode'):
                mode=list(child.text.strip())[0]
                print(mode)
            elif (child.tag=='time'):
                time=int(child.text.strip())
                print(time)
            elif (child.tag=='unsafestate'):
                unsafeState=int(child.text.strip())
                print(unsafeState)
            elif (child.tag=='unsafeparam'):
                unsafeParam=float(child.text.strip())
                print(unsafeParam)
            elif (child.tag=='faults'):
                strFault=child.text.strip()
            elif (child.tag=='A'):
                strA=child.text.strip()
            elif (child.tag=='B'):
                strB=child.text.strip()
            elif (child.tag=='variance'):
                strVar=child.text.strip()
            elif (child.tag=='initialset'):
                strIS=child.text.strip()

        faults=[]
        for f in strFault.split('\n'):
            t=f.strip().split(',')
            faults.append((int(t[0].strip()),int(t[1].strip()),int(t[2].strip()),int(t[3].strip())))
        print(faults)

        variance=[]
        for v in strVar.split('\n'):
            t=v.strip().split(',')
            variance.append((int(t[0].strip()),int(t[1].strip())))
        print(variance)

        initialSet=[]
        for i in strIS.split('\n'):
            t=i.strip().split(',')
            initialSet.append((int(t[0].strip()),int(t[1].strip())))
        print(initialSet)

        a=[]
        for row in strA.split('\n'):
            tmp=row.strip().split(',')
            l=[]
            for t in tmp:
                if (len(t.strip())>0):
                    l.append(int(t.strip()))
            a.append(l)
        A=np.array(a)
        print(A)

        b=[]
        print(strB)
        for row in strB.split('\n'):
            tmp=row.strip().split(',')
            l=[]
            for t in tmp:
                if (len(t.strip())>0):
                    l.append(int(t.strip()))
            b.append(l)
        B=np.array(b)
        print(B)

        #print(int(time/h)+1)
        ConsolidatedVerificationReport.conVerTemplate(A,B,h,mode,faults,variance,time,initialSet,unsafeState,unsafeParam,strName)

    def readFromFileTimeVarying(fname):
        '''A=
        B=
        h=
        mode=
        faults=
        variance=
        initialSet=
        strName=
        time=
        unsafeState=
        unsafeParam=
        Verify.verifyTemplate(A,B,h,mode,faults,variance,initialSet,strName,time,unsafeState,unsafeParam)'''
        tree = ET.parse(fname)
        root = tree.getroot()
        for child in root:
            #print(child.tag,child.text)
            if (child.tag=='h'):
                h=float(child.text.strip())
                print(h)
            elif (child.tag=='name'):
                strName=child.text.strip()
                print(strName)
            elif (child.tag=='mode'):
                mode=list(child.text.strip())[0]
                print(mode)
            elif (child.tag=='time'):
                time=int(child.text.strip())
                print(time)
            elif (child.tag=='unsafestate'):
                unsafeState=int(child.text.strip())
                print(unsafeState)
            elif (child.tag=='unsafeparam'):
                unsafeParam=float(child.text.strip())
                print(unsafeParam)
            elif (child.tag=='faults'):
                strFault=child.text.strip()
            elif (child.tag=='A'):
                strA=child.text.strip()
            elif (child.tag=='B'):
                strB=child.text.strip()
            elif (child.tag=='variance'):
                strVar=child.text.strip()
            elif (child.tag=='initialset'):
                strIS=child.text.strip()
            elif (child.tag=='step'):
                step=int(child.text.strip())

        faults=[]
        for f in strFault.split('\n'):
            t=f.strip().split(',')
            faults.append((int(t[0].strip()),int(t[1].strip()),int(t[2].strip()),int(t[3].strip())))
        print(faults)

        variance=[]
        for v in strVar.split('\n'):
            l=[]
            for t in v.strip().split(';'):
                tmp=t.strip().split(',')
                l.append((int(tmp[0].strip()),int(tmp[1].strip())))
            variance.append(l)
        print(variance)

        initialSet=[]
        for i in strIS.split('\n'):
            t=i.strip().split(',')
            initialSet.append((int(t[0].strip()),int(t[1].strip())))
        print(initialSet)

        a=[]
        for row in strA.split('\n'):
            tmp=row.strip().split(',')
            l=[]
            for t in tmp:
                if (len(t.strip())>0):
                    l.append(int(t.strip()))
            a.append(l)
        A=np.array(a)
        print(A)

        b=[]
        print(strB)
        for row in strB.split('\n'):
            tmp=row.strip().split(',')
            l=[]
            for t in tmp:
                if (len(t.strip())>0):
                    l.append(int(t.strip()))
            b.append(l)
        B=np.array(b)
        print(B)

        #print(int(time/h)+1)
        ConsolidatedVerificationReport.conVerTemplateTimeVarying(A,B,h,mode,faults,variance,time,initialSet,unsafeState,unsafeParam,strName,step)

    def test():
        faults=[(3,3,0,3)]
        #faults=[(2,3,2,3)]
        variance=[(5,5)]
        initialSet=[(-1,1),(-1,1),(-1,1),(-1,1),(-1,1),(-1,1)]
        strName="Example Model"
        time=20
        unsafeState=0
        unsafeParam=0
        ConsolidatedVerificationReport.conVerTemplate(Benchmarks.Test.A,Benchmarks.Test.B,Benchmarks.Test.h,Benchmarks.Test.mode,faults,variance,time,initialSet,unsafeState,unsafeParam,strName)

    def testTimevarying():
        faults=[(3,3,0,3)]
        #faults=[(2,3,2,3)]
        variance=[[(5,5),(1,1)]]
        initialSet=[(-1,1),(-1,1),(-1,1),(-1,1),(-1,1),(-1,1)]
        strName="Example Model (Time-Varying)"
        time=20
        t=1000
        unsafeState=0
        unsafeParam=0
        ConsolidatedVerificationReport.conVerTemplateTimeVarying(Benchmarks.Test.A,Benchmarks.Test.B,Benchmarks.Test.h,Benchmarks.Test.mode,faults,variance,time,initialSet,unsafeState,unsafeParam,strName,t)

    def PKPD2():
        faults=[(0,4,4,1)]
        #faults=[(3,1,4,1)]
        variance=[(5,5)]
        #initialSet=[(1,6),(0,10),(0,10),(1,8),(0,200)]
        initialSet=[(1,6),(0,10),(0,10),(1,8),(1,1)]
        strName="PK/PD Model - II"
        time=20
        unsafeState=0
        unsafeParam=0
        ConsolidatedVerificationReport.conVerTemplate(Benchmarks.PKPD2.A,Benchmarks.PKPD2.B,Benchmarks.PKPD2.h,Benchmarks.PKPD2.mode,faults,variance,time,initialSet,unsafeState,unsafeParam,strName)

    def PKPD2TimeVarying():
        faults=[(0,4,4,1)]
        #faults=[(3,1,4,1)]
        #variance=[[(-1000,1002),(1,1),(0.1,1.9),(-1000,1002),(-1000,1002),(-1,2),(-2,3),(-100,102),(0.8,1.2),(0.5,1.5),(-1,2),(0.1,1.9),(-1,2),(0.5,1.5),(0.6,0.4),(0.8,1.2),(0.7,1.3),(-1,2),(0.8,1.2),(0.5,1.5)]]
        variance=[[(-1000,2000),(-1000,2000),(-1000,2000),(-1000,2000),(-1000,2000),(-1000,2000),(-1000,2000),(-1000,2000),(-1000,2000),(-1000,2000),(-1000,2000),(-1000,2000),(-1000,2000),(-1000,2000),(-1000,2000),(-1000,2000),(-1000,2000),(-1000,2000),(-1000,2000),(-1000,2000),(-1000,2000)]]
        initialSet=[(1,6),(0,10),(0,10),(1,8),(1,1)]
        var=30
        strName="PK/PD Model - II (Time-Varying Perturbation)"
        time=20
        unsafeState=0
        unsafeParam=0
        ConsolidatedVerificationReport.conVerTemplateTimeVarying(Benchmarks.PKPD2.A,Benchmarks.PKPD2.B,Benchmarks.PKPD2.h,Benchmarks.PKPD2.mode,faults,variance,time,initialSet,unsafeState,unsafeParam,strName,var)

    def flightEnvelope():
        faults=[(0,6,6,10)]
        #faults=[(0,8,8,8)]
        variance=[(0.87,1.13)]
        #initialSet=[(3,4),(3,4),(3,5),(3,4),(3,4),(3,5),(-0.50,0.50),(-0.50,0.50),(-0.50,0.50),(-1,1),(-1,1),(-1,1),(5,10),(-0.1,0.1),(-0.1,0.1),(-0.1,0.1)]
        #initialSet=[(3,4),(3,4),(3,5),(3,4),(3,4),(3,5),(-0.50,0.50),(-0.50,0.50),(-0.50,0.50),(-1,1),(-1,1),(-1,1),(1,1),(1,1),(1,1),(1,1)]
        #initialSet=[(3,4),(3,4),(3,5),(3,4),(3,4),(3,5),(-0.50,0.50),(-0.50,0.50),(-0.50,0.50),(-1,1),(1,1),(1,1),(1,1),(1,1),(1,1),(1,1)]
        initialSet=[(3,4),(3,4),(3,4),(3,4),(3,4),(3,5),(-0.5,0.5),(1,1),(1,1),(1,1),(1,1),(1,1),(1,1),(1,1),(1,1),(1,1)]
        strName="Flight Envelope Model"
        time=20
        unsafeState=3
        unsafeParam=-15
        ConsolidatedVerificationReport.conVerTemplate(Benchmarks.FlightEnvelope.A,Benchmarks.FlightEnvelope.B,Benchmarks.FlightEnvelope.h,Benchmarks.FlightEnvelope.mode,faults,variance,time,initialSet,unsafeState,unsafeParam,strName)

    def flightEnvelopeTimeVarying():
        faults=[(0,6,6,10)]
        #faults=[(0,8,8,8)]
        #faults=[(0,2,1,14)]
        #variance=[[(0.99,1.01),(0.98,1.02),(1,1),(0.9,1.1),(0.95,1.05),(0.8,1.2),(0.7,1.3),(0.9,1.1),(1,1),(0.99,1.01),(0.95,1.05),(1,1),(0.98,1.02),(0.8,1.2),(1,1),(0.75,1.25),(1,1),(0.9,1.1),(0.7,1.3),(0.9,1.1)]]
        variance=[[(0.87,1.13),(0.87,1.13),(0.87,1.13),(0.87,1.13),(0.87,1.13),(0.87,1.13),(0.87,1.13),(0.87,1.13),(0.87,1.13),(0.87,1.13),(0.87,1.13),(0.87,1.13),(0.87,1.13),(0.87,1.13),(0.87,1.13),(0.87,1.13),(0.87,1.13),(0.87,1.13),(0.87,1.13),(0.87,1.13)]]
        initialSet=[(3,4),(3,4),(3,4),(3,4),(3,4),(3,5),(-0.50,0.50),(1,1),(1,1),(1,1),(1,1),(1,1),(1,1),(1,1),(1,1),(1,1)]
        strName="Flight Envelope Model (Time-Varying Perturbation)"
        time=20
        var=6
        unsafeState=3
        unsafeParam=-15

        ConsolidatedVerificationReport.conVerTemplateTimeVarying(Benchmarks.FlightEnvelope.A,Benchmarks.FlightEnvelope.B,Benchmarks.FlightEnvelope.h,Benchmarks.FlightEnvelope.mode,faults,variance,time,initialSet,unsafeState,unsafeParam,strName,var)

    def coOPVehiclesI():
        faults=[(0,9,9,1)]
        variance=[(5,5)]
        #initialSet=[(0,0),(0,0),(0,0),(0,0),(0,0),(0,0),(0,0),(0,0),(0,0),(1,1)]
        initialSet=[(1,1),(1,1),(1,1),(1,1),(1,1),(1,1),(1,1),(1,1),(1,1),(1,1)]
        strName="Networked Cooperative Platoon of Vehicles - I"
        time=20
        unsafeState=0
        unsafeParam=-26
        ConsolidatedVerificationReport.conVerTemplate(Benchmarks.CoOPVehiclesI.A,Benchmarks.CoOPVehiclesI.B,Benchmarks.CoOPVehiclesI.h,Benchmarks.CoOPVehiclesI.mode,faults,variance,time,initialSet,unsafeState,unsafeParam,strName)

    def coOPVehiclesII():
        faults=[(0,9,9,1)]
        variance=[(5,5)]
        initialSet=[(1,1),(1,1),(1,1),(1,1),(1,1),(1,1),(1,1),(1,1),(1,1),(1,1)]
        strName="Networked Cooperative Platoon of Vehicles - II"
        time=20
        unsafeState=0
        unsafeParam=0
        ConsolidatedVerificationReport.conVerTemplate(Benchmarks.CoOPVehiclesII.A,Benchmarks.CoOPVehiclesII.B,Benchmarks.CoOPVehiclesII.h,Benchmarks.CoOPVehiclesII.mode,faults,variance,time,initialSet,unsafeState,unsafeParam,strName)

    def coOPVehiclesIITimeVarying():
        faults=[(0,9,9,1)]
        #variance=[[(0.5,1.5),(0.6,1.4),(0.5,1.5),(0.7,1.3),(-1,2),(0.5,1.5),(0.4,0.6),(0.5,1.5),(0.6,1.4),(0.7,1.3),(0.65,1.35),(0.6,1.4),(0.7,1.3),(0.6,1.4),(0.5,1.5),(0.7,1.3),(-1,2),(0.7,1.3),(0.5,1.5),(0.8,1.2)]]
        variance=[[(0.1,1.9),(0.1,1.9),(0.1,1.9),(0.1,1.9),(0.1,1.9),(0.1,1.9),(0.1,1.9),(0.1,1.9),(0.1,1.9),(0.1,1.9),(0.1,1.9),(0.1,1.9),(0.1,1.9),(0.1,1.9),(0.1,1.9),(0.1,1.9),(0.1,1.9),(0.1,1.9),(0.1,1.9),(0.1,1.9)]]
        initialSet=[(1,1),(1,1),(1,1),(1,1),(1,1),(1,1),(1,1),(1,1),(1,1),(1,1)]
        strName="Networked Cooperative Platoon of Vehicles - II (Time-Varying Perturbation)"
        time=20
        var=20
        unsafeState=0
        unsafeParam=1
        ConsolidatedVerificationReport.conVerTemplateTimeVarying(Benchmarks.CoOPVehiclesII.A,Benchmarks.CoOPVehiclesII.B,Benchmarks.CoOPVehiclesII.h,Benchmarks.CoOPVehiclesII.mode,faults,variance,time,initialSet,unsafeState,unsafeParam,strName,var)

    def holes():
        faults=[(0, 4, 4, 4)]
        variance = [(5,5)]
        #initialSet=[(0,0),(0,0),(0,0),(0,0),(0,0),(-100,100),(0,0),(0,0),(-100,100),(-100,100)]
        initialSet=[(0,0),(0,0),(0,0),(0,0),(0,0),(0,0),(-100,-100),(-100,-100)]
        strName="Holes Plant Xp"
        time=20
        unsafeState=0
        unsafeParam=0
        ConsolidatedVerificationReport.conVerTemplate(Benchmarks.HolesPXp.A,Benchmarks.HolesPXp.B,Benchmarks.HolesPXp.h,Benchmarks.HolesPXp.mode,faults,variance,time,initialSet,unsafeState,unsafeParam,strName)

    def motorTransmission1():
        faults=[(0,4,4,3)]
        variance = [(5,5)]
        #initialSet=[(0,0),(-0.08,0.08),(-0.165,-0.165),(-0.01,0.01),(0,0),(70,70),(1,1)]
        initialSet=[(0,0),(-0.08,0.08),(-0.0165,-0.0165),(-0.01,0.01),(0,0),(70,70),(1,1)]
        strName="Motor Transmission - I"
        time=20
        unsafeState=0
        unsafeParam=0
        ConsolidatedVerificationReport.conVerTemplate(Benchmarks.MotorTransmission1.A,Benchmarks.MotorTransmission1.B,Benchmarks.MotorTransmission1.h,Benchmarks.MotorTransmission1.mode,faults,variance,time,initialSet,unsafeState,unsafeParam,strName)

    def motorTransmission1TimeVarying():
        faults=[(0,4,4,3)]
        #variance = [[(0.7,1.3),(-0.2,1),(0.01,1.99),(-0.1,1),(0.8,1.2),(-0.3,1),(0.7,1.3),(1,1),(0.8,1.2),(0.6,0.4),(0.65,0.35),(0.8,1.2),(1,1),(1,1.2),(0.7,1),(0.8,1),(1,1.1),(1,1.2),(1,1),(0.9,1.1)]]
        variance = [[(-0.01,1),(-0.01,1),(-0.01,1),(-0.01,1),(-0.01,1),(-0.01,1),(-0.01,1),(-0.01,1),(-0.01,1),(-0.01,1),(-0.01,1),(-0.01,1),(-0.01,1),(-0.01,1),(-0.01,1),(-0.01,1),(-0.01,1),(-0.01,1),(-0.01,1),(-0.01,1)]]
        #initialSet=[(0,0),(-0.08,0.08),(-0.165,-0.165),(-0.01,0.01),(0,0),(70,70),(1,1)]
        initialSet=[(0,0),(-0.08,0.08),(-0.0165,-0.0165),(-0.01,0.01),(0,0),(70,70),(1,1)]
        var=2
        strName="Motor Transmission - I (Time-Varying Perturbation)"
        time=20
        unsafeState=2
        unsafeParam=-0.02
        ConsolidatedVerificationReport.conVerTemplateTimeVarying(Benchmarks.MotorTransmission1.A,Benchmarks.MotorTransmission1.B,Benchmarks.MotorTransmission1.h,Benchmarks.MotorTransmission1.mode,faults,variance,time,initialSet,unsafeState,unsafeParam,strName,var)

    def motorTransmission2():
        faults=[(2,3,0,2)]
        variance = [(5,5)]
        initialSet=[(0,0),(-0.08,0.08),(-0.165,-0.165),(-0.01,0.01),(0,0)]
        strName="Motor Transmission - II"
        time=20
        unsafeState=0
        unsafeParam=0
        ConsolidatedVerificationReport.conVerTemplate(Benchmarks.MotorTransmission2.A,Benchmarks.MotorTransmission2.B,Benchmarks.MotorTransmission2.h,Benchmarks.MotorTransmission2.mode,faults,variance,time,initialSet,unsafeState,unsafeParam,strName)
