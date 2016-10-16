def frontforward(A,B,PI,O,t):  
    # 0<=t<T时，输出向量，t=T时输出最终概率
    # n是状态数，m是观测数
    n,m = np.shape(B)
    T = np.shape(O)[0]
    if 0< t < T: # 当t小于等于0时，按t等于0计算
        front = frontforward(A,B,PI,O,t-1)
        return [np.dot(front,A[:,i])*B[i][O[t]] for i in range(n)]
    elif t <= 0:
        return [PI[i]*B[i][O[t]] for i in range(n)]
    else: #当t大于等于T时，按t等于length计算
        return np.sum(frontforward(A,B,PI,O,t-1))
        
        
        
def backforward(A,B,PI,O,t):
    # 0<=t<=T-1时输出向量，t=-1时输出最终概率
    n,m = np.shape(B)
    T = np.shape(O)[0]
    if 0<=t<=T-2:
        back = backforward(A,B,PI,O,t+1)
        return [np.dot(back*B[:,O[t+1]],A[i,:]) for i in range(n)]
    elif t >= T-1:
        return [1]*n
    else:
        return np.sum(backforward(A,B,PI,O,0)*PI*B[:,O[0]])
        
        
def backfront(A,B,PI,O,t):
    n,m = np.shape(B)
    T = np.shape(O)[0]
    if 0<=t<=T-2:
        return np.sum([frontforward(A,B,PI,O,t)[i]*backforward(A,B,PI,O,t+1)[j]*A[i,j]*B[j,O[t+1]] for i in range(n) for j in range(n)])
    elif t >= T-1:
        return frontforward(A,B,PI,O,T)
    else:
        return backforward(A,B,PI,O,-1)
        
        
def gammat(A,B,PI,O,t,i):
    prosingle = frontforward(A,B,PI,O,t)[i]*backforward(A,B,PI,O,t)[i]
    prosum = np.sum([frontforward(A,B,PI,O,t)[i]*backforward(A,B,PI,O,t)[i] for j in range(n)])
    return 1.0*prosingle/prosum
    
    
def xit(A,B,PI,O,t,i,j):
    prosingle = frontforward(A,B,PI,O,t)[i]*backforward(A,B,PI,O,t+1)[j]*A[i,j]*B[j][O[t+1]]
    prosum = np.sum([frontforward(A,B,PI,O,t)[i]*backforward(A,B,PI,O,t+1)[j]*A[i,j]*B[j][O[t+1]] for i in range(n) for j in range(n)])
    return 1.0*prosingle/prosum
    
    
print '后向算法backforward的测试：'
for t in range(5,-2,-1):
    print 't=%d:'%t,backforward(A,B,PI,O,t)
    
    
print '计算概率backfront的测试：'
for t in range(-1,6):
    print 't=%d:'%t,backfront(A,B,PI,O,t)
