# -*- coding: utf-8 -*-
"""
Created on Thu Jan 17 16:15:59 2019

@author: xiaoy
"""
def Outage_schedule(C,q,r,delta_x):
    '''
    C:容量
    --------------
    FOR:强迫停运率
    --------------
    r:强迫停运时间
    --------------
    delta_x:步长
    --------------
    '''
    mu=1/r
    I=C//delta_x+1
    return [[C-delta_x*i for i in range(I)],\
            [q if i!=0 else 1.0000 for i in range(I)],\
            [q*mu if i!=0 else 0.0000 for i in range(I)]]
def Parallel_Outage_schedule(a,b,delta_x):
    '''
    a:原件a的停运表
    --------------
    b:原件b的停运表
    --------------
    '''
    NI=len(a[0])+len(b[0])-1
    xi=list(range(a[0][0]+b[0][0],a[0][-1]+b[0][-1]-delta_x,-delta_x))
    a_p=[i-j for i,j in zip(a[1],(a[1]+[0])[1::])]
    a_f=[i-j for i,j in zip(a[2],(a[2]+[0])[1::])]
    pi_x=list()
    fi_x=list()
    for k in range(NI):
        pp=0
        pf=0
        for i,pa in enumerate(a_p[::]):
            if len(b[1])-1>=(k-i)>=0:
                pp=b[1][k-i]*pa+pp
                pf=b[2][k-i]*pa+b[1][k-i]*a_f[i]+pf
            elif k-i>0:
                pp=0.0*pa+pp
                pf=0.0*pa+0.0*a_f[i]+pf
            elif k-i<0:
                pp=1.0*pa+pp
                pf=0.0*pa+1.0*a_f[i]+pf
        pi_x.append(pp)
        fi_x.append(pf)
    return [xi,pi_x,fi_x]
def Series_Outage_schedule(aa,bb):
    '''
    aa:原件a的停运表
    --------------
    bb:原件b的停运表
    --------------
    '''
    import copy
    a=copy.deepcopy(aa)
    b=copy.deepcopy(bb)
    while a[0][0]>b[0][0]:
        del a[0][0],a[1][0],a[2][0]
    else:
        a[1][0]=1.0
        a[2][0]=0.0
    while a[0][0]<b[0][0]:
        del b[0][0],b[1][0],b[2][0]
    else:
        b[1][0]=1.0
        b[2][0]=0.0
    xi=copy.deepcopy(a[0])
    pi_x=list()
    fi_x=list()
    for k,l in enumerate(range(len(a[0]))):
        pi_x.append(a[1][k]+b[1][k]-a[1][k]*b[1][k])
        fi_x.append(a[2][k]*(1-b[1][k])+(1-a[1][k])*b[2][k])
    return [xi,pi_x,fi_x]
def Line_Degree_Form(x,Line,Line_Degree,Keypoint):
    ''' 
    形成线路级别表的函数
    '''
    for y in Line[:]:
        if Keypoint in y:
            z=list(y)
            z.remove(Keypoint)
            z=z[0]
            if x in Line_Degree:
                Line_Degree[x].append((Keypoint,z))
            else:
                Line_Degree[x]=[(Keypoint,z)]
            Line.remove(y)
            Line_Degree_Form(x+1,Line,Line_Degree,z)