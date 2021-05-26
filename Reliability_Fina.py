# -*- coding: utf-8 -*-
"""
Created on Thu Jan 17 17:05:28 2019

@author: xiaoy
"""
from Reliability_Func import *
Point_Reliability=dict()
Line_Reliability=dict()
Line_Degree=dict()
with open('Power_Information.txt','r') as f:
    data=f.readlines()
data1=data[0][:].strip()
data1=[int(x) for x in data1.split(' ')]
N_P=data1[0]
N_L=data1[1]
N_D=data1[2]
delta_x=data1[3]
Key_Point=data1[4]
del data[0],data1
for i in range(N_P):
    dataI=data[i][:]
    dataI=[float(x) if (y!=0)and(y!=1) else int(x) for y,x in enumerate(dataI.split(' '))]
    if dataI[0] not in Point_Reliability:        
        Point_Reliability[dataI[0]]=Outage_schedule(dataI[1],dataI[2],dataI[3],delta_x)
    else:
        Point_Reliability[dataI[0]]=Parallel_Outage_schedule(Outage_schedule(dataI[1],dataI[2],dataI[3],delta_x),\
                         Point_Reliability[dataI[0]],delta_x)
del data[:N_P],dataI
for i in range(N_L):
    dataI=data[i][:]
    dataI=[float(x) if (y!=0)and(y!=1)and(y!=2) else int(x) for y,x in enumerate(dataI.split(' '))]
    if ((dataI[0],dataI[1]) not in Line_Reliability) and ((dataI[1],dataI[0]) not in Line_Reliability):        
        Line_Reliability[(dataI[0],dataI[1])]=Outage_schedule(dataI[2],dataI[3],dataI[4],delta_x)
    elif (dataI[0],dataI[1]) in Line_Reliability:
        Line_Reliability[(dataI[0],dataI[1])]=Parallel_Outage_schedule(Outage_schedule(dataI[2],dataI[3],dataI[4],delta_x),\
                         Line_Reliability[(dataI[0],dataI[1])],delta_x)
    elif (dataI[1],dataI[0]) in Line_Reliability:
        Line_Reliability[(dataI[1],dataI[0])]=Parallel_Outage_schedule(Outage_schedule(dataI[2],dataI[3],dataI[4],delta_x),\
                         Line_Reliability[(dataI[1],dataI[0])],delta_x)
del data[:N_L],dataI
for i in range(N_D):
    I=int(data[0][5])
    del data[0]
    dataI=[[],[],[]]
    while(data[0][0:5]!='point'):       
        dataII=[float(x) if y!=0 else int(x) for y,x in enumerate(data[0].split(' '))]
        del data[0]
        dataI[0].append(dataII[0])
        dataI[1].append(dataII[1])
        dataI[2].append(dataII[2])
        if len(data)==0:
            break
    if I not in Point_Reliability: 
        Point_Reliability[I]=dataI
    else:
        Point_Reliability[I]=Parallel_Outage_schedule(dataI,\
                 Point_Reliability[I],delta_x)
del data,dataI,dataII,I
Line=list(Line_Reliability)
Line_Degree_Form(0,Line,Line_Degree,Key_Point)
del Line
Line_Reversed=list(Line_Degree)
Line_Reversed.sort(reverse=True)
for i in Line_Reversed:
    for j in Line_Degree[i]:
        if j[0] not in Point_Reliability:
            if j in Line_Reliability:
                Point_Reliability[j[0]]=Series_Outage_schedule(Line_Reliability[j],Point_Reliability[j[1]])
            else:
                Point_Reliability[j[0]]=Series_Outage_schedule(Line_Reliability[(j[1],j[0])],Point_Reliability[j[1]])
        else:
            if j in Line_Reliability:
                Point_Reliability[j[0]]=Parallel_Outage_schedule(Series_Outage_schedule(Line_Reliability[j],Point_Reliability[j[1]]),\
                                 Point_Reliability[j[0]],delta_x)
            else:
                Point_Reliability[j[0]]=Parallel_Outage_schedule(Series_Outage_schedule(Line_Reliability[(j[1],j[0])],Point_Reliability[j[1]]),\
                                 Point_Reliability[j[0]],delta_x)
del Line_Reversed,Line_Degree,Line_Reliability,i,j
P_D=[Point_Reliability[Key_Point][1][k] for k,x in enumerate(Point_Reliability[Key_Point][0]) if x<0]
LOLP=P_D[0]
EENS=sum(P_D)*delta_x*24
LOLF=[Point_Reliability[Key_Point][2][k] for k,x in enumerate(Point_Reliability[Key_Point][0]) if x<0][0]
LOLE=LOLP*24
D=LOLE/LOLF