def coor_get_lammps(file,*positional_parameters, **keyword_parameters):
    xyz=open(file)
    natom=0
    nbond=0
    nangle=0
    atype=0
    flag_use=1
    split_line=[]
    flags_out=0

    if ("images" in keyword_parameters): flags_out=1


    if ("noflag" in keyword_parameters): flag_use=0

    while natom==0:
        read_line=xyz.readline()
        split_line=read_line.split()
#        print split_line
        if len(split_line)==2 and split_line[1]=="atoms":
            natom=int(split_line[0])

    nbond_found=0
    while nbond_found==0 and ('nbond' in keyword_parameters):	
        read_line=xyz.readline()
        split_line=read_line.split()
#        print split_line
        if len(split_line)==2 and split_line[1]=="bonds":
            nbond=int(split_line[0])
            nbond_found=1

    nangle_found=0
    while nangle_found==0 and ('nangle' in keyword_parameters):
        read_line=xyz.readline()
        split_line=read_line.split()
#        print split_line
        if len(split_line)==2 and split_line[1]=="angles":
            nangle=int(split_line[0])
            nangle_found=1

    natype_found=0
    while natype_found==0 and ('masses' in keyword_parameters):
        read_line=xyz.readline()
        split_line=read_line.split()
#        print split_line
        if len(split_line)==3 and split_line[1]=="atom":	
            atype=int(split_line[0])
            natype_found=1
    vel_yes=0
    if ('velocities' in keyword_parameters):
#        print("SSSSUUUUUUUPPPPPP") 
        while vel_yes!=1:
            read_line=xyz.readline()
            split_line=read_line.split()
            if len(split_line)>1:
                if split_line[0]=="Velocities": vel_yes=1
        output=[]
        velocities=[[] for i in range(4)]
        split_line=[]
        while len(split_line)==0:
            read_line=xyz.readline()
            split_line=read_line.split()
        while len(split_line)>0:    
            for i in range(3): velocities[i].append(float(split_line[i+1]))
            velocities[3].append(float(split_line[0]))
            read_line=xyz.readline()
            split_line=read_line.split()     
        
        output.append(velocities)        
        return  output


    tilts_found=0
    tilts=[]
    break_count=0
    while tilts_found==0 and ('tilts' in keyword_parameters):

#        print "finding tilts: "
        read_line=xyz.readline()
        split_line=read_line.split()
#        if len(split_line)==6: 
#            print "tilts finding", split_line
#            print split_line[3]
        if len(split_line)==6 and split_line[3]=="xy":
            tilts=[0,0,0]
            print("split= ", split_line)
            for i in range(3): tilts[i]=float(split_line[i]) 
            tilts_found=1
        else:
          break_count+=1
        if break_count==100: break


#    print "tilts= ",tilts

    offset=0
    if ('q' in keyword_parameters):
#        print("getting charges")
        offset=1

    mol=0
    if ("mol" in keyword_parameters):
        mol=1


#    print("offset= ", offset)
#    print("we got there")
#    print("natom= ",natom)



    split_line=[]
    output=[]
    if ('natom' in keyword_parameters):
#        print("natom keyword found")
        output.append(natom)
        if ('nbond' in keyword_parameters):
            output.append(nbond)
        if ('nangle' in keyword_parameters):
            output.append(nangle)	

    elif ('tilts' in keyword_parameters):
#        print("getting tilts")
        output.append(tilts)    
        

    


    elif ('masses' in keyword_parameters):
#        print("getting masses")
        while split_line!=["Masses"]:
            read_line=xyz.readline()
            split_line=read_line.split()
#            print split_line
        read_line=xyz.readline()
        masses=[0 for i in range(atype)]
#        print("atype= ",atype)
        for i in range(atype):
            read_line=xyz.readline()
            split_line=read_line.split()
            masses[int(split_line[0])-1]=float(split_line[1])

        output.append(masses)





    else:
        xlen=0
        while xlen==0:
            read_line=xyz.readline()
            split_line=read_line.split()
            if len(split_line)==4:
                if split_line[2]=='xlo': 
                    xlen=float(split_line[1])-float(split_line[0])
        read_line=xyz.readline()
        split_line=read_line.split()
        ylen=float(split_line[1])-float(split_line[0])
        read_line=xyz.readline()
        split_line=read_line.split()
        zlen=float(split_line[1])-float(split_line[0])


        print("xlen, ylen, zlen: ", xlen, ylen, zlen)

        Afound=0
        while Afound==0:
            read_line=xyz.readline()
            split_line=read_line.split()
            if len(split_line)>0:
                if split_line[0]=="Atoms": Afound=1
        read_line=xyz.readline()
        x=[0 for i in range(natom)]
        y=[0 for i in range(natom)]
        z=[0 for i in range(natom)]
        mols=[0 for i in range(natom)]
        type=[0 for i in range(natom)]
        charges=[0 for i in range(natom)]
        xflag=[0 for i in range(natom)]
        yflag=[0 for i in range(natom)]
        zflag=[0 for i in range(natom)]
#        print("get coors")
        for i in range(natom):
            read_line=xyz.readline()
            split_line=read_line.split()
            id=int(split_line[0])
#            xflag=0
#            yflag=0
#            zflag=0
            if len(split_line)==9+offset and flag_use==1:
                xflag[id-1]=int(split_line[6+offset])
                yflag[id-1]=int(split_line[7+offset])
                zflag[id-1]=int(split_line[8+offset])
   
#            print("split_line= ", split_line) 
            x[id-1]=float(split_line[3+offset])+xflag[id-1]*xlen
            y[id-1]=float(split_line[4+offset])+yflag[id-1]*ylen
            z[id-1]=float(split_line[5+offset])+zflag[id-1]*zlen
            mols[id-1]=int(split_line[1])
#            print("new coors= ", x[id-1], y[id-1], z[id-1])
#            print split_line
#            print x[id-1], y[id-1], z[id-1]
 
            type[id-1]=int(split_line[2])
            charges[id-1]=float(split_line[3])

        if flags_out==1:
            output=[xflag,yflag,zflag]

        elif offset==1:
            output=[x,y,z,type,mols,charges,[xlen,ylen,zlen]]
#            print(output)
        else:
            output=[x,y,z,type,mols,[xlen,ylen,zlen]]
#    print(output)
    xyz.close()
    return output







def coor_get_lammps_image(file,*positional_parameters, **keyword_parameters):
    xyz=open(file)
    natom=0
    nbond=0
    nangle=0
    atype=0
    flag_use=0
    split_line=[]
    if ("noflag" in keyword_parameters): flag_use=0

    while natom==0:
        read_line=xyz.readline()
        split_line=read_line.split()
#        print split_line
        if len(split_line)==2 and split_line[1]=="atoms":
            natom=int(split_line[0])

    nbond_found=0
    while nbond_found==0 and ('nbond' in keyword_parameters):
        read_line=xyz.readline()
        split_line=read_line.split()
#        print split_line
        if len(split_line)==2 and split_line[1]=="bonds":
            nbond=int(split_line[0])
            nbond_found=1

    nangle_found=0
    while nangle_found==0 and ('nangle' in keyword_parameters):
        read_line=xyz.readline()
        split_line=read_line.split()
#        print split_line
        if len(split_line)==2 and split_line[1]=="angles":
            nangle=int(split_line[0])
            nangle_found=1

    natype_found=0
    while natype_found==0 and ('masses' in keyword_parameters):
        read_line=xyz.readline()
        split_line=read_line.split()
#        print split_line
        if len(split_line)==3 and split_line[1]=="atom":
            atype=int(split_line[0])
            natype_found=1
    vel_yes=0
    if ('velocities' in keyword_parameters):
#        print("SSSSUUUUUUUPPPPPP") 
        while vel_yes!=1:
            read_line=xyz.readline()
            split_line=read_line.split()
            if len(split_line)>1:
                if split_line[0]=="Velocities": vel_yes=1
        output=[]
        velocities=[[] for i in range(4)]
        split_line=[]
        while len(split_line)==0:
            read_line=xyz.readline()
            split_line=read_line.split()
        while len(split_line)>0:
            for i in range(3): velocities[i].append(float(split_line[i+1]))
            velocities[3].append(float(split_line[0]))
            read_line=xyz.readline()
            split_line=read_line.split()

        output.append(velocities)
        return  output


    tilts_found=0
    tilts=[0,0,0]
    while tilts_found==0 and ('tilts' in keyword_parameters):
#        print "finding tilts: "
        read_line=xyz.readline()
        split_line=read_line.split()
#        if len(split_line)==6: 
#            print "tilts finding", split_line
#            print split_line[3]
        if len(split_line)==6 and split_line[3]=="xy":
            for i in range(3): tilts[i]=float(split_line[i])
            tilts_found=1

#    print "tilts= ",tilts

    offset=0
    if ('q' in keyword_parameters):
#        print("getting charges")
        offset=1

    mol=0
    if ("mol" in keyword_parameters):
        mol=1


#    print("offset= ", offset)
#    print("we got there")
#    print("natom= ",natom)



    split_line=[]
    output=[]
    if ('natom' in keyword_parameters):
#        print("natom keyword found")
        output.append(natom)
        if ('nbond' in keyword_parameters):
            output.append(nbond)
        if ('nangle' in keyword_parameters):
            output.append(nangle)

    elif ('tilts' in keyword_parameters):
#        print("getting tilts")
        output.append(tilts)





    elif ('masses' in keyword_parameters):
#        print("getting masses")
        while split_line!=["Masses"]:
            read_line=xyz.readline()
            split_line=read_line.split()
#            print split_line
        read_line=xyz.readline()
        masses=[0 for i in range(atype)]
#        print("atype= ",atype)
        for i in range(atype):
            read_line=xyz.readline()
            split_line=read_line.split()
            masses[int(split_line[0])-1]=float(split_line[1])

        output.append(masses)



    else:
        xlen=0
        while xlen==0:
            read_line=xyz.readline()
            split_line=read_line.split()
            if len(split_line)==4:
                if split_line[2]=='xlo':
                    xlen=float(split_line[1])-float(split_line[0])
        read_line=xyz.readline()
        split_line=read_line.split()
        ylen=float(split_line[1])-float(split_line[0])
        read_line=xyz.readline()
        split_line=read_line.split()
        zlen=float(split_line[1])-float(split_line[0])


        print("xlen, ylen, zlen: ", xlen, ylen, zlen)

        Afound=0
        while Afound==0:
            read_line=xyz.readline()
            split_line=read_line.split()
            if len(split_line)>0:
                if split_line[0]=="Atoms": Afound=1
        read_line=xyz.readline()
        x=[0 for i in range(natom)]
        y=[0 for i in range(natom)]
        z=[0 for i in range(natom)]
        mols=[0 for i in range(natom)]
        type=[0 for i in range(natom)]
        charges=[0 for i in range(natom)]
#        print("get coors")
        for i in range(natom):
            read_line=xyz.readline()
            split_line=read_line.split()
            id=int(split_line[0])
            xflag=0
            yflag=0
            zflag=0
            if len(split_line)==9+offset and flag_use==1:
                xflag=int(split_line[6+offset])
                yflag=int(split_line[7+offset])
                zflag=int(split_line[8+offset])

#            print("split_line= ", split_line)
            x[id-1]=float(split_line[3+offset])+xflag*xlen
            y[id-1]=float(split_line[4+offset])+yflag*ylen
            z[id-1]=float(split_line[5+offset])+zflag*zlen
            mols[id-1]=int(split_line[1])
#            print("new coors= ", x[id-1], y[id-1], z[id-1])
#            print split_line
#            print x[id-1], y[id-1], z[id-1]

            type[id-1]=int(split_line[2])
            charges[id-1]=float(split_line[3])

        if offset==1:
            output=[x,y,z,type,mols,charges,[xlen,ylen,zlen]]
#            print(output)
        else:
            output=[x,y,z,type,mols,[xlen,ylen,zlen]]
#    print(output)
    xyz.close()
    return output












def coor_get_lammps_mol(file,*positional_parameters, **keyword_parameters):
    xyz=open(file)
    natom=0
    nbond=0
    nangle=0
    atype=0
    while natom==0:

        read_line=xyz.readline()
        split_line=read_line.split()
        if len(split_line)==2 and split_line[1]=="atoms":
            natom=int(split_line[0])



    while nbond==0 and ('nbond' in keyword_parameters):     
        read_line=xyz.readline()
        split_line=read_line.split()
        if len(split_line)==2 and split_line[1]=="bonds":
            nbond=int(split_line[0])


    while nangle==0 and ('nangle' in keyword_parameters):
        read_line=xyz.readline()
        split_line=read_line.split()
        if len(split_line)==2 and split_line[1]=="angles":
            nangle=int(split_line[0])

    while atype==0 and ('masses' in keyword_parameters):
        read_line=xyz.readline()
        split_line=read_line.split()
        if len(split_line)==3 and split_line[1]=="atom":
            atype=int(split_line[0])





    offset=0
    if ('q' in keyword_parameters):
        if keyword_parameters['q']=="True":
            
            offset=1
    print("we got there")
    print("natom= ",natom)



    split_line=[]
    output=[]
    if ('natom' in keyword_parameters):
        print("natom keyword found")
        output.append(natom)
        if ('nbond' in keyword_parameters):
            output.append(nbond)
        if ('nangle' in keyword_parameters):
            output.append(nangle)


    elif ('masses' in keyword_parameters):
        print("getting masses")
        while split_line!=["Masses"]:
            read_line=xyz.readline()
            split_line=read_line.split()
        read_line=xyz.readline()
        masses=[0 for i in range(atype)]
        for i in range(atype):
            read_line=xyz.readline()
            split_line=read_line.split()
            masses[int(split_line[0])-1]=float(split_line[1])

        output.append(masses)

                        
    else:
        print("suppppp")
        xlen=0
        while xlen==0:
            read_line=xyz.readline()
            split_line=read_line.split()
        #               print split_line
            if len(split_line)==4:
                print(split_line)
                if split_line[2]=='xlo': xlen=float(split_line[1])-float(split_line[0])
        #                       print xlen
        #               if xlen!=0: break

        #       xlen=float(split_line[1])-float(split_line[3])
                read_line=xyz.readline()
                split_line=read_line.split()
                ylen=float(split_line[1])-float(split_line[0])
                read_line=xyz.readline()
                split_line=read_line.split()
                zlen=float(split_line[1])-float(split_line[0])

        Afound=0
        print("Finding Atoms")

        while Afound==0:
            read_line=xyz.readline()
            split_line=read_line.split()
            if len(split_line)>0:
                if split_line[0]=="Atoms": 
                    Afound=1
                    print("Afound")

        read_line=xyz.readline()
        x=[0 for i in range(natom)]
        y=[0 for i in range(natom)]
        z=[0 for i in range(natom)]
        type=[0 for i in range(natom)]
        mol=[0 for i in range(natom)]
        charges=[0 for i in range(natom)]
        for i in range(natom):
            read_line=xyz.readline()
            split_line=read_line.split()
            id=int(split_line[0])

            x[id-1]=float(split_line[3+offset])
            y[id-1]=float(split_line[4+offset])
            z[id-1]=float(split_line[5+offset])
            type[id-1]=int(split_line[2])
            mol[id-1]=int(split_line[1])
            charges[id-1]=float(split_line[3])

        if offset==1:
            output=[x,y,z,type,charges,mol,[xlen,ylen,zlen]]
        else:
            output=[x,y,z,type,mol,[xlen,ylen,zlen]]

    xyz.close()
    return output









def coor_get_lammps_unordered(file,*positional_parameters, **keyword_parameters):
    xyz=open(file)
    natom=0
    nbond=0
    nangle=0
    atype=0
    flag_use=1

    if ("noflag" in keyword_parameters): flag_use=0


    while natom==0:
        read_line=xyz.readline()
        split_line=read_line.split()
#        print split_line
        if len(split_line)==2 and split_line[1]=="atoms":
            natom=int(split_line[0])

    while nbond==0 and ('nbond' in keyword_parameters):
        read_line=xyz.readline()
        split_line=read_line.split()
#        print split_line
        if len(split_line)==2 and split_line[1]=="bonds":
            nbond=int(split_line[0])


    while nangle==0 and ('nangle' in keyword_parameters):
        read_line=xyz.readline()
        split_line=read_line.split()
#        print split_line
        if len(split_line)==2 and split_line[1]=="angles":
            nangle=int(split_line[0])

    while atype==0 and ('masses' in keyword_parameters):
        read_line=xyz.readline()
        split_line=read_line.split()
#        print split_line
        if len(split_line)==3 and split_line[1]=="atom":
            atype=int(split_line[0])

    tilts_found=0
    tilts=[0,0,0]
    while tilts_found==0 and ('tilts' in keyword_parameters):
#        print "finding tilts: "
        read_line=xyz.readline()
        split_line=read_line.split()
#        if len(split_line)==6: 
#            print "tilts finding", split_line
#            print split_line[3]
        if len(split_line)==6 and split_line[3]=="xy":
            for i in range(3): tilts[i]=float(split_line[i])
            tilts_found=1

#    print "tilts= ",tilts



    vel_yes=0
    if ('velocities' in keyword_parameters):
#        print("SSSSUUUUUUUPPPPPP")
        while vel_yes!=1:
            read_line=xyz.readline()
            split_line=read_line.split()
#            print(split_line)
            if len(split_line)>0:
                if split_line[0]=="Velocities": vel_yes=1
#        print("Sup2")
        output=[]
        velocities=[[] for i in range(4)]
        split_line=[]
        while len(split_line)==0:
            read_line=xyz.readline()
            split_line=read_line.split()
        while len(split_line)>0:
            for i in range(3): velocities[i].append(float(split_line[i+1]))
            velocities[3].append(int(split_line[0]))
            read_line=xyz.readline()
            split_line=read_line.split()

        output.append(velocities)
        return  output









    offset=0
    if ('q' in keyword_parameters):
        print("getting charges")
        offset=1

    mol=0
    if ("mol" in keyword_parameters):
        mol=1

    print("offset= ", offset)
    print("we got there")
    print("natom= ",natom)



    split_line=[]
    output=[]
    if ('natom' in keyword_parameters):
        print("natom keyword found")
        output.append(natom)
        if ('nbond' in keyword_parameters):
            output.append(nbond)
        if ('nangle' in keyword_parameters):
            output.append(nangle)

    elif ('tilts' in keyword_parameters):
        print("getting tilts")
        output.append(tilts)





    elif ('masses' in keyword_parameters):
        print("getting masses")
        while split_line!=["Masses"]:
            read_line=xyz.readline()
            split_line=read_line.split()
#            print split_line
        read_line=xyz.readline()
        masses=[0 for i in range(atype)]
        print("atype= ",atype)
        for i in range(atype):
            read_line=xyz.readline()
            split_line=read_line.split()
            masses[int(split_line[0])-1]=float(split_line[1])

        output.append(masses)





    else:
        xlen=0
        while xlen==0:
            read_line=xyz.readline()
            split_line=read_line.split()
            if len(split_line)==4:
                if split_line[2]=='xlo':
                    xlen=float(split_line[1])-float(split_line[0])
        read_line=xyz.readline()
        split_line=read_line.split()
        ylen=float(split_line[1])-float(split_line[0])
        read_line=xyz.readline()
        split_line=read_line.split()
        zlen=float(split_line[1])-float(split_line[0])


        print("xlen, ylen, zlen: ", xlen, ylen, zlen)

        Afound=0
        while Afound==0:
            read_line=xyz.readline()
            split_line=read_line.split()
            if len(split_line)>0:
                if split_line[0]=="Atoms": Afound=1
        read_line=xyz.readline()
        split_line=read_line.split()
        while len(split_line)==0:
            read_line=xyz.readline()
            split_line=read_line.split()

        x=[0 for i in range(natom)]
        y=[0 for i in range(natom)]
        z=[0 for i in range(natom)]
        mols=[0 for i in range(natom)]
        type=[0 for i in range(natom)]
        charges=[0 for i in range(natom)]
        print("get coors")
        id=[0 for i in range(natom)]

        for i in range(natom):
            id[i]=int(split_line[0])
            xflag=0
            yflag=0
            zflag=0
            if len(split_line)==9+offset and flag_use==1:
                xflag=int(split_line[6+offset])
                yflag=int(split_line[7+offset])
                zflag=int(split_line[8+offset])

            x[i]=float(split_line[3+offset])+xflag*xlen
            y[i]=float(split_line[4+offset])+yflag*ylen
            z[i]=float(split_line[5+offset])+zflag*zlen
            mols[i]=int(split_line[1])

            type[i]=int(split_line[2])
            charges[i]=float(split_line[3])

            read_line=xyz.readline()
            split_line=read_line.split()

        if offset==1:
            output=[x,y,z,type,mols,charges,id,[xlen,ylen,zlen]]
#            print(output)
        else:
            output=[x,y,z,type,mols,id,[xlen,ylen,zlen]]
#    print(output)
    xyz.close()
    return output
