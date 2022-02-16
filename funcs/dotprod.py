def	dotprod(v1,v2):



	l1=sum([v1[i]*v1[i] for i in range(len(v1))])**.5
	l2=sum([v2[i]*v2[i] for i in range(len(v1))])**.5

	

	dot_sum=sum([v1[i]*v2[i] for i in range(len(v1))])
	dot=dot_sum/l1/l2

	return	dot
	
