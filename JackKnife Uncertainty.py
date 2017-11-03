from __future__ import division
from math import*
import numpy as np

endProgram=False
while endProgram==False:
    values=[]
    endLoop=False
    while endLoop==False:
        val=float(raw_input("Enter value: "))
        values.append(val)
        option=raw_input("End input? Enter yes or no.")
        if option=="yes":
            endLoop=True

    AvgVal=sum(values)/len(values)
    summ=0
    for j in range(0,len(values)):
        summ=summ+(values[j]-AvgVal)**2
    print sqrt(((len(values)-1)/len(values))*summ)
    option2=raw_input("End Program? Enter yes or no.")
    if option2=="yes":
        endProgram=True
    
