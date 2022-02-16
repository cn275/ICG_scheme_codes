def	gen_rot(coors,a,b,d,*positional_parameters,**keyword_parameters):
	import math,lammpstrjconvert
	new_coors=[[0 for i in range(len(coors[0]))] for j in range(3)]
	for i in range(len(coors[0])):
		for j in range(3):
			coors[j][i]=float(coors[j][i])

	ca=math.cos(a)
	cb=math.cos(b)
	cd=math.cos(d)
	
	sa=math.sin(a)
	sb=math.sin(b)
	sd=math.sin(d)

	for i in range(len(coors[0])):
		new_coors[0][i]=(coors[0][i])*ca*cb+coors[1][i]*(ca*sb*sd-sa*cd)+coors[2][i]*(ca*sb*cd+sa*sd)
		new_coors[1][i]=(coors[0][i])*sa*cb+coors[1][i]*(sa*sb*sd+ca*cd)+coors[2][i]*(sa*sb*cd-sd*ca)
		new_coors[2][i]=-(coors[0][i])*sb+coors[1][i]*(cb*sd)+coors[2][i]*(cd*cb)

	
	if "compare" in keyword_parameters:
		if keyword_parameters['compare']=="True":
			comptrj=[[i for i in j] for j in coors]
			for i in range(3): comptrj[i].extend(new_coors[i])
			types=[1+i/len(coors[0]) for i in range(2*len(coors[0]))]	
			lammpstrjconvert.lammpstrjconvert(comptrj[0],comptrj[1],comptrj[2],[i+1 for i in range(len(comptrj[0]))],types,1,1,1,"comp_rot.lammpstrj")

	return	new_coors	
