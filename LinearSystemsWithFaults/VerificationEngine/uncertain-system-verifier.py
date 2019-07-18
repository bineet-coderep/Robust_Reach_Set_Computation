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

Please refer to the README.md to use or modify this Verification Engine.

Link to the repository: <>
'''

from RobustReachSetAPI import *


# Please call APIs as per how-to-use.md here:

#Verify.flightEnvelope()


#-------------------------------------------------------




# Please comment out this section if you are calling APIs

if (len(sys.argv)>=2):
    ConsolidatedVerificationReport.readFromFile(sys.argv[1])
else:
    print("Nothing to execute, please call a function!")

#--------------------------------------------
