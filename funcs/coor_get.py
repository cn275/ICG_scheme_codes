def coor_get(file,key, *positional_parameters, **keyword_parameters):
    xyz=open(file)
    split_line=[]
#    print(file, key)

    unwrap=0
    unscaled=1
    offset=0
#    print(keyword_parameters)
    if ("idoffset" in keyword_parameters):
        print("offset= ",int(keyword_parameters["idoffset"])) 
        offset=int(keyword_parameters["idoffset"])
    if ("unwrap" in keyword_parameters):
        unwrap=int(keyword_parameters["unwrap"])
    if ("unscaled" in keyword_parameters):
        unscaled=int(keyword_parameters["unscaled"])


    for i in range(2):
        read_line=xyz.readline()
        split_line=read_line.split()	
#        print(split_line)
    time=int(split_line[0])
#    print("time= ",time)

    for i in range(2):
        read_line=xyz.readline()
#        print(read_line)
        split_line=read_line.split()
#    print("natm= ",read_line)
    natom=int(read_line.split()[0])
    told=int(time)


#    print(int(key)==int(time))
#    print(int(key)==time)

    if int(time) !=int(key):	
        while split_line != [str(key)]:
            read_line=xyz.readline()
            split_line=read_line.split()


    while len(split_line) <6:
        read_line=xyz.readline()
        split_line=read_line.split()
    read_line=xyz.readline()
    split_line=read_line.split()
#    print(split_line)


    xy=0
    xz=0
    yz=0	
    xlen=(float(split_line[1])-float(split_line[0]))
    if ('tilts' in keyword_parameters): xy=float(split_line[2])
    read_line=xyz.readline()
    split_line=read_line.split()
    if ('tilts' in keyword_parameters): xz=float(split_line[2])
    ylen=(float(split_line[1])-float(split_line[0]))
    read_line=xyz.readline()
    split_line=read_line.split()
    if ('tilts' in keyword_parameters): yz=float(split_line[2])
    zlen=(float(split_line[1])-float(split_line[0]))
    x=[0 for i in range(natom)]
    y=[0 for i in range(natom)]
    z=[0 for i in range(natom)]
    id=[0 for i in range(natom)]
    type=[0 for i in range(natom)]
    read_line=xyz.readline()
 
    for i in range(natom):
        read_line=xyz.readline()
        split_line=read_line.split()
#        print(read_line)
        id=int(split_line[0])-offset
        x[id-1]=float(split_line[2-unwrap])*float(xlen)**unscaled
        y[id-1]=float(split_line[3-unwrap])*float(ylen)**unscaled
        z[id-1]=float(split_line[4-unwrap])*float(zlen)**unscaled
        type[id-1]=int(float(split_line[1]))

    if unwrap==1:
        type=[0 for i in type]
    xyz.close()
    return [x,y,z,type,[xlen,ylen,zlen,xy,xz,yz]]







def coor_get_last(file,natom, *positional_parameters, **keyword_parameters):


    xyz=open(file)
    emptycount=0
    key=0

    last_found=0
    read_line=xyz.readline()
    read_line=xyz.readline()
    split_line=read_line.split()
    key=int(split_line[0])
    while last_found==0:
        for i in range(natom+9):
            read_line=xyz.readline()
        split_line=read_line.split()
        if len(split_line)>0:
            key=int(split_line[0])
#            print "key= ", key
        else:
            last_found=1

    xyz.close()


#    print "last_frame: ", key





    xyz=open(file)
    split_line=[]
#    print(file, key)

    unwrap=0
    unscaled=1
    offset=0
#    print(keyword_parameters)
    if ("idoffset" in keyword_parameters):
#        print("offset= ",int(keyword_parameters["idoffset"]))
        offset=int(keyword_parameters["idoffset"])
    if ("unwrap" in keyword_parameters):
        unwrap=int(keyword_parameters["unwrap"])
    if ("unscaled" in keyword_parameters):
        unscaled=int(keyword_parameters["unscaled"])


    for i in range(2):
        read_line=xyz.readline()
        split_line=read_line.split()
#        print(split_line)
    time=int(split_line[0])
#    print("time= ",time)


    if int(time) !=int(key):
        while split_line != [str(key)]:
            for i in range(natom+9):
                read_line=xyz.readline()
            split_line=read_line.split()
#            print split_line 

    for i in range(2):
        read_line=xyz.readline()
#        print(read_line)
        split_line=read_line.split()
#    print("natm= ",read_line)
    natom=int(read_line.split()[0])
    told=int(time)


    while len(split_line) <6:
        read_line=xyz.readline()
        split_line=read_line.split()
    read_line=xyz.readline()
    split_line=read_line.split()
#    print(split_line)


    xy=0
    xz=0
    yz=0
    xlen=(float(split_line[1])-float(split_line[0]))
    if ('tilts' in keyword_parameters): xy=float(split_line[2])
    read_line=xyz.readline()
    split_line=read_line.split()
    if ('tilts' in keyword_parameters): xz=float(split_line[2])
    ylen=(float(split_line[1])-float(split_line[0]))
    read_line=xyz.readline()
    split_line=read_line.split()
    if ('tilts' in keyword_parameters): yz=float(split_line[2])
    zlen=(float(split_line[1])-float(split_line[0]))
    x=[0 for i in range(natom)]
    y=[0 for i in range(natom)]
    z=[0 for i in range(natom)]
    id=[0 for i in range(natom)]
    type=[0 for i in range(natom)]


    read_line=xyz.readline()

    for i in range(natom):
        read_line=xyz.readline()
        split_line=read_line.split()
#        print(read_line)
        id=int(split_line[0])-offset
        x[id-1]=float(split_line[2-unwrap])*float(xlen)**unscaled
        y[id-1]=float(split_line[3-unwrap])*float(ylen)**unscaled
        z[id-1]=float(split_line[4-unwrap])*float(zlen)**unscaled
        type[id-1]=int(float(split_line[1]))

    if unwrap==1:
        type=[0 for i in type]
    xyz.close()
    return [x,y,z,type,[xlen,ylen,zlen,xy,xz,yz]]





















def coor_get_nozero(file,key):
    xyz=open(file)
    split_line=[]
    while split_line != [str(key)]:
        read_line=xyz.readline()
        split_line=read_line.split()
    read_line=xyz.readline()
    read_line=xyz.readline()
    natom=int(read_line)
    while len(split_line) <6:
        read_line=xyz.readline()
        split_line=read_line.split()


    read_line=xyz.readline()
    split_line=read_line.split()
    print(split_line)
    xlen=(float(split_line[1])-float(split_line[0]))
    read_line=xyz.readline()
    split_line=read_line.split()
    ylen=(float(split_line[1])-float(split_line[0]))
    read_line=xyz.readline()
    split_line=read_line.split()
    zlen=(float(split_line[1])-float(split_line[0]))
    x=[0 for i in range(natom)]
    y=[0 for i in range(natom)]
    z=[0 for i in range(natom)]
    id=[0 for i in range(natom)]
    type=[0 for i in range(natom)]
    read_line=xyz.readline()
    empty=0
    real=0

    for i in range(natom):
        read_line=xyz.readline()
        split_line=read_line.split()
        id=int(split_line[0])
        if (id==0) or (float(split_line[1]))==0:
            empty+=1	
        else:
            real+=1
            x[id-1]=float(split_line[2])*float(xlen)
            y[id-1]=float(split_line[3])*float(ylen)
            z[id-1]=float(split_line[4])*float(zlen)
            type[id-1]=int(split_line[1])

    empty=type.count(0)
    print("empty= ",empty)
    print("real= ",real)
    x=[x[i] for i in range(natom-empty)]
    y=[y[i] for i in range(natom-empty)]
    z=[z[i] for i in range(natom-empty)]
    type=[type[i] for i in range(natom-empty)]
    xyz.close()
    print("returned coor")
    return [x,y,z,type,[xlen,ylen,zlen]]





#def coor_get_aa(file,lfile,key,bondmax,nchain,chainl):
#	import bond_get
#	import pbcorrectsingle
#	import coor_get_lammps
#        xyz=open(file)
#	bonds=bond_get.bond_get(lfile)
#        split_line=[]
#        while split_line != [str(key)]:
#                read_line=xyz.readline()
#                split_line=read_line.split()
#        read_line=xyz.readline()
#        read_line=xyz.readline()
#        natom_real=coor_get_lammps.coor_get_lammps(lfile,natom="yes")
#	natom=int(read_line)
#        while split_line != ["ITEM:", "BOX","BOUNDS", "pp", "pp", "pp"]:
#                read_line=xyz.readline()
#                split_line=read_line.split()
#        read_line=xyz.readline()
#        split_line=read_line.split()
#        print split_line
#        xlen=(float(split_line[1])-float(split_line[0]))
#        read_line=xyz.readline()
#        split_line=read_line.split()
#        ylen=(float(split_line[1])-float(split_line[0]))
#        read_line=xyz.readline()
#        split_line=read_line.split()
#        zlen=(float(split_line[1])-float(split_line[0]))
#	id_list=[]
#	x=[]
#	y=[]
#	z=[]
#	type=[]
#        read_line=xyz.readline()
#        empty=0
#        real=0
#        for i in range(natom):
#                read_line=xyz.readline()
#                split_line=read_line.split()
#                id=int(split_line[0])
#		if id>0:
#			x.append(float(split_line[2])*float(xlen))
#			y.append(float(split_line[3])*float(ylen))
#			z.append(float(split_line[4])*float(zlen))
#			type.append(int(split_line[1]))
#			id_list.append(id)
#
#	xreal=[0 for i in range(natom_real)]
#        yreal=[0 for i in xreal]
#        zreal=[0 for i in xreal]
#        typereal=[69 for i in xreal]
#        idreal=[0 for i in xreal]
#
#	for i in range(nchain):
#		chainid=[1+j+i*chainl for j in range(chainl)]
#		ncop=[id_list.count(j) for j in chainid]
#		print chainid, ncop
#		ids=[[k for k in range(len(id_list)) if id_list[k]==j] for j in chainid]		
#		print ids
#		print " "
#		print " "
#		chainfounds=[0 for i in range(chainl)]
#		for j in range(chainl):
#			if len(ids[j])==1: 
#				print "One coordinate option for ",chainid[j], x[ids[j][0]],y[ids[j][0]],z[ids[j][0]]
#				xreal[chainid[j]-1]=x[ids[j][0]]
#				yreal[chainid[j]-1]=y[ids[j][0]]
#				zreal[chainid[j]-1]=z[ids[j][0]]
#				typereal[chainid[j]-1]=type[ids[j][0]]
#				chainfounds[j]=1
#
#		while sum(chainfounds)<chainl:
#			for j in range(chainl):
#				if chainfounds[j]==0:
#					print "searching for real coordinates of %d"%chainid[j]
#					if j==0 and chainfounds[1]==1:
#						print "chain founds partner= ",chainid[1], xreal[chainid[1]-1],yreal[chainid[1]-1],zreal[chainid[1]-1]
#						chainfounds[0]=1
#						mindist=2*xlen
#						minid=0
#						for k in range(len(ids[j])):
#							print "possible coordinate= ",x[ids[j][k]],y[ids[j][k]],z[ids[j][k]]
#
#							dist=pbcorrectsingle.pbcorrect(x[ids[j][k]],y[ids[j][k]],z[ids[j][k]],xreal[chainid[1]-1],yreal[chainid[1]-1],zreal[chainid[1]-1],xlen,ylen,zlen)[-1]				
#							if dist<mindist:
#								mindist=dist
#								minid=ids[j][k]
#						xreal[chainid[0]-1]=x[minid]
#						yreal[chainid[0]-1]=y[minid]
#						zreal[chainid[0]-1]=z[minid]
#						typereal[chainid[0]-1]=type[minid]
#						print "found coor= ",xreal[chainid[0]-1], yreal[chainid[0]-1], zreal[chainid[0]-1]
#
#					if j+1==chainl and chainfounds[chainl-2]==1:
#						mindist=2*xlen
#                                                minid=0
#						print "chain founds partner= ",chainid[chainl-2], xreal[chainid[chainl-2]-1],yreal[chainid[chainl-2]-1],zreal[chainid[chainl-2]-1]
#						chainfounds[chainl-1]=1				
#						for k in range(len(ids[j])):
#                                                        print "possible coordinate= ",x[ids[j][k]],y[ids[j][k]],z[ids[j][k]]
#
#                                                        dist=pbcorrectsingle.pbcorrect(x[ids[j][k]],y[ids[j][k]],z[ids[j][k]],xreal[chainid[chainl-2]-1],yreal[chainid[chainl-2]-1],zreal[chainid[chainl-2]-1],xlen,ylen,zlen)[-1]
#                                                        if dist<mindist:
#                                                                mindist=dist
#                                                                minid=ids[j][k]
#                                                xreal[chainid[j]-1]=x[minid]
#                                                yreal[chainid[j]-1]=y[minid]
#                                                zreal[chainid[j]-1]=z[minid]
#						typereal[chainid[0]-1]=type[minid]
#						print 	"found coor= ", xreal[chainid[j]-1], yreal[chainid[j]-1], zreal[chainid[j]-1]
#
#
#					if j>0 and j+1<chainl:
#						mindist=2*xlen
#                                                minid=0				
#						if chainfounds[j+1]==1:
#							print "chain founds partner= ",chainid[j+1], xreal[chainid[j+1]-1],yreal[chainid[j+1]-1],zreal[chainid[j+1]-1]
#							chainfounds[j]=1
#							for k in range(len(ids[j])):
#	                                                        print "possible coordinate= ",x[ids[j][k]],y[ids[j][k]],z[ids[j][k]]
#	                                                        dist=pbcorrectsingle.pbcorrect(x[ids[j][k]],y[ids[j][k]],z[ids[j][k]],xreal[chainid[j+1]-1],yreal[chainid[j+1]-1],zreal[chainid[j+1]-1],xlen,ylen,zlen)[-1]
#                                                        	if dist<mindist:
#                                                                	mindist=dist
#                                                                	minid=ids[j][k]
#                                                	xreal[chainid[j]-1]=x[minid]
#                                                	yreal[chainid[j]-1]=y[minid]
#                                                	zreal[chainid[j]-1]=z[minid]
#							typereal[chainid[0]-1]=type[minid]
#							print   "found coor= ", xreal[chainid[j]-1], yreal[chainid[j]-1], zreal[chainid[j]-1]
#
#
#						if chainfounds[j-1]==1:
#							chainfounds[j]=1
#							print "chain founds partner= ",chainid[j-1], xreal[chainid[j-1]-1],yreal[chainid[j-1]-1],zreal[chainid[j-1]-1]
#							for k in range(len(ids[j])):
#                                                                print "possible coordinate= ",x[ids[j][k]],y[ids[j][k]],z[ids[j][k]]
#                                                                dist=pbcorrectsingle.pbcorrect(x[ids[j][k]],y[ids[j][k]],z[ids[j][k]],xreal[chainid[j-1]-1],yreal[chainid[j-1]-1],zreal[chainid[j-1]-1],xlen,ylen,zlen)[-1]
#                                                        	if dist<mindist:
#                                                        	        mindist=dist
#                                                        	        minid=ids[j][k]
#                                                	xreal[chainid[j]-1]=x[minid]
#                                                	yreal[chainid[j]-1]=y[minid]
#                                                	zreal[chainid[j]-1]=z[minid]		
#							typereal[chainid[0]-1]=type[minid]
#							print   "found coor= ", xreal[chainid[j]-1], yreal[chainid[j]-1], zreal[chainid[j]-1]
#
#
#        xyz.close()
#	output=[xreal,yreal,zreal,typereal,[xlen,ylen,zlen]]
#        return output

