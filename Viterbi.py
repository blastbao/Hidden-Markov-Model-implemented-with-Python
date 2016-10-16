def viterbi(A,B,PI,O):
    T = np.shape(O)[0]
    N,K = np.shape(B)  
    I = np.array([np.nan]*T) 
    delta = np.zeros((T,N))
    phi = np.zeros((T,N))
    for i in range(N):
        delta[0,i] = logSpace.elnproduct(logSpace.eln(PI[i]),logSpace.eln(B[i,O[0]]))
        phi[0,i] = 0
    for t in range(1,T):
        for i in range(N):
            maxvalue = logSpace.elnproduct(delta[t-1,0],logSpace.eln(A[0,i]))
            maxindex = 0
            for j in range(1,N):
                candidvalue = logSpace.elnproduct(delta[t-1,j],logSpace.eln(A[j,i]))
                if maxvalue < candidvalue:
                    maxvalue = candidvalue
                    maxindex = j
            delta[t,i] = logSpace.elnproduct(maxvalue,logSpace.eln(B[i,O[t]]))
            phi[t,i] = maxindex
    P = np.max(delta[-1,:])
    I[-1] = np.where(delta[-1,:] == P)[0][0]
    P = logSpace.eexp(P)
    for t in range(T-2,-1,-1):
        I[t] = phi[t+1][I[t+1]]
    return P,I
    
    
    
A = np.array([[0.5,0.2,0.3],
              [0.3,0.5,0.2],
              [0.2,0.3,0.5]])
B = np.array([[0.5,0.5],
              [0.4,0.6],
              [0.7,0.3]])
PI = np.array([0.2,0.4,0.4])
O = np.array([0,1,0])
T = np.shape(O)[0]
N,K = np.shape(B)


P,I = viterbi(A,B,PI,O)
print '最优路径：'
print I
print '该最优路径的概率：'
print P
