def avg_std(A):
        N=len(A)
        sum1=0.0
        for i in range(N):
                sum1=sum1+float(A[i])
        Average=sum1/N
        sum2=0.0
        for i in range(N):
                sum2=sum2+(float(A[i])-Average)**2
        SD=(sum2/(N-1))**.5
        return Average, SD

