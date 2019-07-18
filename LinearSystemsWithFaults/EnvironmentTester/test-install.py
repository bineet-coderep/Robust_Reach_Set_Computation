'''
Environment Test:

Run the following:
python test-install.py

If everything is installed correctly, then this
should run without any error and display
"Environment is Ready!!" at the end.
'''

import numpy as np
import itertools
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

        #print(U)
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


        model.addConstr(rVars[self.unsafeVariable]<=self.unsafeParam)

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

    class Test:
        A=np.array([
        [1,0],
        [0,0],
        ])
        B=np.array([
        ])
        h=0.01
        mode='.'

class Verify:

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
        dyn=Dynamics(inp.formMatA(),n,f)
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

class VerifyFaultFree:

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

class ConsolidatedVerificationReport:
    #Under Construction
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

    def installTest():
        faults=[(0,1,1,1)]
        variance=[(5,5)]
        initialSet=[(-1,1),(-1,1)]
        strName="Install Tester Model"
        time=20
        unsafeState=0
        unsafeParam=0
        ConsolidatedVerificationReport.conVerTemplate(Benchmarks.Test.A,Benchmarks.Test.B,Benchmarks.Test.h,Benchmarks.Test.mode,faults,variance,time,initialSet,unsafeState,unsafeParam,strName)


ConsolidatedVerificationReport.installTest()
print("\n\n=======================")
print("Environment is Ready!!")
print("=======================\n\n")
