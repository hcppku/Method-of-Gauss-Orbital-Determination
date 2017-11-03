#excercise 2 hw2

#part a:
from __future__ import division
from math import*
import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt

#returns sum of every element in array
def N(array):
    sum=0
    for rowIndex in range(0,len(array)):
        for columnIndex in range(0,len(array[0])):
            sum=sum+array[rowIndex,columnIndex]
    return sum

#returns array of sums of each column
def sumColumn(array):
    sumColumn=[]
    for columnIndex in range(0,len(array[0])):
        sum=0
        for rowIndex in range(0,len(array)):
            sum=sum+array[rowIndex,columnIndex]
        sumColumn.append(sum)
    return sumColumn
    
#returns array of sums of each row
def sumRow(array):
    sumRow=[]
    for rowIndex in range(0,len(array)):
        sum=0
        for columnIndex in range(0,len(array[0])):
            sum=sum+array[rowIndex,columnIndex]
        sumRow.append(sum)
    return sumRow

#returns x position for centroid
def xPosition(array):
    sc = sumColumn(array)
    sumBeforeMean=0
    divideSum=0
    
    for scIndex in range(0,len(sc)):
        sumBeforeMean=sumBeforeMean+scIndex*sc[scIndex]
        divideSum=divideSum+sc[scIndex]

    return sumBeforeMean/divideSum

#returns y position for centroid
def yPosition(array):
    sr=sumRow(array)
    sumBeforeMean=0
    divideSum=0
    for srIndex in range(0,len(sr)):
        sumBeforeMean=sumBeforeMean+srIndex*sr[srIndex]
        divideSum=divideSum+sr[srIndex]
    return sumBeforeMean/divideSum

#returns x and y coordinate for centroid
def centroid (array):
    return xPosition(array),yPosition(array)

#returns array of distance squared of each column to the x centroid coordinate
def distX(array):
    distanceX=[]
    xPos=xPosition(array)
    for column in range(0,len(array[0])):
        distanceX.append((xPos-column)**2)
    return distanceX

#returns array of distance squared of each column to the x centroid coordinate
def distY(array):
    distanceY=[]
    yPos=yPosition(array)
    for row in range(0,len(array)):
        distanceY.append((yPos-row)**2)
    return distanceY

#finds stdev for x coordinate
def stdevX(array):
    dX=distX(array)
    sCol=sumColumn(array)
    summation=0.0
    for col in range(0,len(array[0])):
        summation=summation+dX[col]*sCol[col]
    N1=N(array)/1000
    N2=(N(array)-1)/1000
    result=(1/1000)*sqrt(summation/(N1*N2))
    return result

#finds stdev for y coordinate
def stdevY(array):
    dY=distY(array)
    sRow=sumRow(array)
    summation=0
    for row in range(0,len(array)):
        summation=summation+dY[row]*sRow[row]
    N1=N(array)/1000
    N2=(N(array)-1)/1000
    result=(1/1000)*sqrt(summation/(N1*N2))
    return result

#test:
#testArray=np.array([[0,33,21,33,8],[0,56,51,53,26],[23,120,149,73,18],[55,101,116,50,16],[11,78,26,2,10]])
#print("centroid coordinates (x,y): ")
#print(centroid(testArray))
#print("\nstdev of x: ")
#print(stdevX(testArray))
#print("\nstdev of y: ")
#print(stdevY(testArray))


#part b:
#image=fits.getdata("Obs1(59).fit")

#displays image
#def display(image):
    #plt.imshow(image)
    #plt.gray()
    #plt.show()

#return centroid, uncertainties of image
def centroidIm(image, width, height, x, y):
    cutImage=image[y-int(height/2):y+int(height/2)+1,x-int(width/2):x+int(width/2)+1]
    uncX=stdevX(cutImage)
    uncY=stdevY(cutImage)
    cent=centroid(cutImage)
    newCent=[cent[0]+x-int(width/2),cent[1]+y-int(height/2)]
    return newCent,uncX,uncY

#test:
##test=centroidIm(image,7,7,358,477)
##print("\n\ncentroid of asteroid: ")
##print(test[0])
##print("\nuncertainty x: ")
##print(test[1])
##print("\nuncertainty y: ")
##print(test[2])
##
##test=centroidIm(image,5,5,358,477)
##print("\n\ncentroid of asteroid: ")
##print(test[0])
##print("\nuncertainty x: ")
##print(test[1])
##print("\nuncertainty y: ")
##print(test[2])
#display(image)

loop=True
while loop==True:
    fileName=fits.getdata(raw_input("Enter image file name: "))
    xy=float(raw_input("Enter xy dimensions of centroid box: "))
    x=float(raw_input("Enter approximate x: "))
    y=float(raw_input("Enter approximate y: "))
    print xy
    cent=centroidIm(fileName,xy,xy,x,y)
    print("\n\ncentroid of asteroid: ")
    print(cent[0])
    print("\nuncertainty x: ")
    print(cent[1])
    print("\nuncertainty y: ")
    print(cent[2])

    choice=int(raw_input("Enter 1 for another calculation or 2 to end: "))
    if choice==2:
        loop=False

    

    
