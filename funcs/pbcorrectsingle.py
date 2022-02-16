def pbcorrect(x1,y1,z1,x2,y2,z2,xlen,ylen,zlen):
	xvec=x1-x2
	yvec=y1-y2
	zvec=z1-z2	
	rorig=(xvec**2+yvec**2+zvec**2)**.5 
	if xvec>.5*xlen:
		xvec=xvec-float(xlen)
	if xvec<-.5*xlen:
		xvec=xvec+float(xlen)
	if yvec>.5*ylen:
		yvec=yvec-float(ylen)
	if yvec<-.5*ylen:
		yvec=yvec+float(ylen)
	if zvec>.5*zlen:
		zvec=zvec-float(zlen)
	if zvec<-.5*zlen:
		zvec=zvec+float(zlen)
	return [xvec,yvec,zvec,(xvec**2+yvec**2+zvec**2)**.5]


def get_flag(x1,y1,z1,x2,y2,z2,xlen,ylen,zlen):

    flags=[0,0,0]
    xvec=x1-x2
    yvec=y1-y2
    zvec=z1-z2
    rorig=(xvec**2+yvec**2+zvec**2)**.5

    if xvec>.5*xlen:
         while xvec>.5*xlen:
            xvec=xvec-float(xlen)
            flags[0]+=1
    
    if xvec<-.5*xlen:
        while xvec>.5*xlen:
            xvec=xvec+float(xlen)
            flags[0]-=1

    if yvec>.5*ylen:
        while yvec>.5*ylen:
            yvec=yvec-float(ylen)
            flags[1]+=1
    if yvec<-.5*ylen:
        while yvec<-.5*ylen:
            yvec=yvec+float(ylen)
            flags[1]-=1

    if zvec>.5*zlen:
        while zvec>.5*zlen:
            zvec=zvec-float(zlen)
            flags[2]+=1
    if zvec<-.5*zlen:
        while zvec<-.5*zlen:
            zvec=zvec+float(zlen)
            flags[2]-=1
    return [xvec,yvec,zvec,flags,(xvec**2+yvec**2+zvec**2)**.5]
