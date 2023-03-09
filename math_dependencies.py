import numpy as np
from scipy.linalg import circulant
from OpenGL.GL import *
import glfw
phi=(1+np.sqrt(5))/2
def getcoords(expr):
    #expr is some tuple that looks like (0,0,1,1,2)
    expr=np.array(expr,dtype=int)
    i = np.argmax(expr)
    if i == 4:
        return -expr[1:4]+1-expr[0]
    else:
        c=circulant([+(phi-1),+phi,1])/2
        t=expr[np.array([i^3,i^1,i^2],dtype=int)]-1+expr[4]
        temp = -c@t
        h=((np.nonzero(t)[0][0])-1%3)^i^2
        j=1-2*np.array([(h&2)//2,h&1,(h&1)^((h&2)//2)])
        return temp*j
def get_transform(perm):
    basis=[np.array([0,0,1,1,2]),np.array([0,1,0,1,2]),np.array([0,1,1,0,2])]
    return np.array([getcoords(b[perm]) for b in basis]).T

def get_transform_4d(perm):
    A=np.zeros((4,4))
    A[:3,:3]=get_transform(perm)
    A[3][3]=1
    return A
import itertools
perms3d = dict(zip(itertools.permutations([1,2,3]),range(6)))


def get_cycle(cycle):
    p1=np.roll(cycle,-np.where(cycle==0)[0][0]%5)
    t=-int(np.where(p1==4)[0][0])%5
    return perms3d[(p1[t],p1[2*t%5],p1[3*t%5])]
    
def get_tri_colours(perm):
    l=[
        np.array([0,2,1,3,4]),
        np.array([0,1,3,2,4]),
        np.array([0,3,2,1,4])]
    
    return [get_cycle(perm[l[i]]) for i in range(3)]

def get_pent_colours(perm):
    l=[
        np.array([0,2,1,3,4]),
        np.array([1,0,4,2,3]),
        np.array([4,1,3,0,2]),
        np.array([3,4,2,1,0]),
        np.array([2,3,0,4,1]),]
    return [get_cycle(perm[l[i]]) for i in range(5)]
def sgn(perm):
    s=perm.shape[0]
    A=np.zeros((s,s))
    A[range(s),perm]=1
    return np.linalg.det(A)

def lp3(i,j):
    ans=np.zeros(5,dtype=int)
    ans[0]=i
    ans[4]=j
    ans[1:4]=np.nonzero((np.arange(5)!=i)*(np.arange(5)!=j))[0]
    if sgn(ans)<0:
        ans[2],ans[3]=ans[3],ans[2]
    return ans

def lp5(i,j):
    ans=np.zeros(5,dtype=int)
    ans[0]=i
    ans[1]=j
    ans[4]=4
    ans[2:4]=np.nonzero((np.arange(4)!=i)*(np.arange(4)!=j))[0]
    if sgn(ans)<0:
        ans[2],ans[3]=ans[3],ans[2]
    return ans

def rotation_matrix(axis,theta):
    axis=axis/np.linalg.norm(axis)
    mat=np.cos(theta)*np.eye(3)+ np.sin(theta)*(np.cross(axis,np.eye(3)))+axis*axis.reshape((-1,1))*(1-np.cos(theta))
    ans=np.zeros((4,4))
    ans[:3,:3]=mat
    ans[3,3]=1
    return ans

tsymlist=[lp3(i%5,(i//5+i+1)%5) for i in range(20)]
tcidlist=np.array([get_tri_colours(s) for s in tsymlist],dtype=int)
ttransformlist=[get_transform_4d(s) for s in tsymlist]

psymlist=[lp5(i%4,(i//4+i+1)%4) for i in range(12)]

pcidlist=np.array([get_pent_colours(s) for s in psymlist],dtype=int)

ptransformlist=[get_transform_4d(s) for s in psymlist]

