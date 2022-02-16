def lammpsinput (x,y,z,id,types,*positional_parameters,**keyword_parameters):

    if ('bonds' in keyword_parameters):	bonds=keyword_parameters['bonds']
    else:					bonds=[[],[],[]]
    nbonds=len(bonds[0])
    if ('angles' in keyword_parameters): 	angles=keyword_parameters['angles']
    else:					angles=[[],[],[],[]]
    nangles=len(angles[0])	
    if ('dihedrals' in keyword_parameters):	dihedrals=keyword_parameters['dihedrals']

    else:					dihedrals=[[],[],[],[],[]]
    ndihedrals=len(dihedrals[0])

    if ('impropers' in keyword_parameters): impropers=keyword_parameters['impropers']
    else:                                   impropers=[[],[],[],[],[]]
    nimpropers=len(impropers[0])

    if ('charges' in keyword_parameters): 	charges=keyword_parameters['charges']
    else:					charges=["" for i in x]


    if ('lens' in keyword_parameters):	lens=keyword_parameters['lens']
    else:					lens=[max(x)-min(x),max(y)-min(y),max(z)-min(z)]
    if ('filename' in keyword_parameters):  filename=keyword_parameters['filename'] 
    else:					filename="Generic_Name.txt"

    if ('extra_lines' in  keyword_parameters): extra_lines=keyword_parameters['extra_lines']
    else:   extra_lines=[]


    if ('tilts' in keyword_parameters):	tilts=keyword_parameters['tilts']
    else:					tilts=[]

    vel_yes=0
    if ("velocities" in keyword_parameters): 
        vel_yes=1
        velocities=keyword_parameters['velocities']

    if ("images" in keyword_parameters): [ximage,yimage,zimage]=keyword_parameters['images']


    natom=len(x)
#    print("natom= ",natom)

    outfile = open(filename,'w')
    outfile.write('LAMMPS Description\n')
    outfile.write('\n')
    outfile.write(str(len(x)))
    outfile.write(' atoms\n')
    outfile.write(str(len(bonds[0])))
    outfile.write(' bonds\n')
    outfile.write(str(len(angles[0])))
    outfile.write(' angles\n')
    outfile.write(str(len(dihedrals[0])))
    outfile.write(' dihedrals\n')
    outfile.write(str(len(impropers[0])))
    outfile.write(' impropers\n')

    outfile.write('\n')



    if ("masses" in keyword_parameters):
        outfile.write('%s atom types\n'%len(keyword_parameters["masses"]))
    else:
        outfile.write('%s atom types\n'%max(types))

#    print "bonds= ", bonds[0]
#    print "angles= ", angles[0]
#    print "dihedrals= ", dihedrals[0]


    if nbonds>0:		
        
        if ("bondcoeff" in keyword_parameters):
            outfile.write('%s bond types\n'%len(keyword_parameters['bondcoeff']))
        else:
            extra_btype_offset=0
            if ("extrabond" in keyword_parameters): extra_btype_offset=int(keyword_parameters['extrabond']) 
            outfile.write('%s bond types\n'%(max([int(i) for i in bonds[0]])+extra_btype_offset))


    if nangles>0:		

        if ("anglecoeff" in keyword_parameters):
            outfile.write('%s angle types\n'%len(keyword_parameters['anglecoeff']))
        else:
            outfile.write('%s angle types\n'%max([int(i) for i in angles[0]]))


    if ndihedrals>0:	

        if ("dihedralcoeff" in keyword_parameters):
            outfile.write('%s dihedral types\n'%len(keyword_parameters['dihedralcoeff']))
        else:
            outfile.write('%s dihedral types\n'%max([int(i) for i in dihedrals[0]]))

    if nimpropers>0:	
    
        if ("impropercoeff" in keyword_parameters):
            outfile.write('%s improper types\n'%len(keyword_parameters['impropercoeff']))
        else:
            outfile.write('%s improper types\n'%max([int(i) for i in impropers[0]]))	
    
    outfile.write('\n')
    for i in extra_lines: outfile.write("%s\n"%i)
    outfile.write('\n')

    outfile.write('\n')
    outfile.write('0')
    outfile.write(' ')
    outfile.write(str(lens[0]))
    outfile.write(' ')
    outfile.write('xlo ')
    outfile.write('xhi\n')
    outfile.write(str(0))
    outfile.write(' ')
    outfile.write(str(lens[1]))
    outfile.write(' ')
    outfile.write('ylo ')
    outfile.write('yhi\n')
    outfile.write(str(0))
    outfile.write(' ')
    outfile.write(str(lens[2]))
    outfile.write(' ')
    outfile.write('zlo ')
    outfile.write('zhi\n')
    outfile.write('\n')

    if ('tilts' in keyword_parameters):
        outfile.write('%s %s %s xy xz yz\n'%(tilts[0],tilts[1],tilts[2]))
        outfile.write('\n')

    outfile.write('Masses\n')
    outfile.write('\n')
	
    if ("masses" in keyword_parameters):
        masses=keyword_parameters['masses']	
        print("masses= ",masses)
        for i in range(len(masses)):
            outfile.write("%s %s\n"%(i+1,masses[i]))	

    else:
        for i in range(max(types)):
            outfile.write('%s 1.0\n'%(i+1))

    outfile.write('\n')


    if (nbonds>0) and ('paircoeff' in keyword_parameters):
        outfile.write('Pair Coeffs\n')
        outfile.write('\n')
        if ('paircoeff' in keyword_parameters):
            bcount=0
            for i in keyword_parameters['paircoeff']:
                bcount+=1
                outfile.write("%s %s\n"%(bcount,i))
        outfile.write('\n')


    if (nbonds>0) and ('bondcoeff' in keyword_parameters):	
        outfile.write('Bond Coeffs\n')
        outfile.write('\n')
        if ('bondcoeff' in keyword_parameters):
            bcount=0
            for i in keyword_parameters['bondcoeff']:
                bcount+=1
                outfile.write("%s %s\n"%(bcount,i))
        outfile.write('\n')

    if (nangles>0) and ('anglecoeff' in keyword_parameters):
        outfile.write('Angle Coeffs\n')
        outfile.write('\n')
        #now we can use loops to do the heavy lifting
        if ('anglecoeff' in keyword_parameters):
            acount=0
            for i in keyword_parameters['anglecoeff']:
                acount+=1
                outfile.write("%s %s\n"%(acount,i))
            outfile.write('\n')

    if (ndihedrals>0) and ('dihedralcoeff' in keyword_parameters):
        outfile.write('\n')
        outfile.write('Dihedral Coeffs\n')
        outfile.write('\n')
        #now we can use loops to do the heavy lifting
        if ('dihedralcoeff' in keyword_parameters):
            dcount=0
            for i in keyword_parameters['dihedralcoeff']:
                dcount+=1
                outfile.write("%s %s\n"%(dcount,i))
        outfile.write('\n')


    if (nimpropers>0) and ('impropercoeff' in keyword_parameters):
        outfile.write('\n')
        outfile.write('Improper Coeffs\n')
        outfile.write('\n')
        #now we can use loops to do the heavy lifting
        if ('impropercoeff' in keyword_parameters):
            icount=0
            for i in keyword_parameters['impropercoeff']:
                icount+=1
                outfile.write("%s %s\n"%(icount,i))
        outfile.write('\n')





    outfile.write("Atoms\n\n")

    print("len types= ",len(types))
	
    if ("molecules" in keyword_parameters): 
        molecules=keyword_parameters["molecules"]
#        print "molecules..... bitches"
    else:					molecules=[1 for i in x]

    print("len(id), len(molecules), len(types), len(charges), len(x), len(y), len(z)")
    print(len(id), len(molecules), len(types), len(charges), len(x), len(y), len(z))

    if ("rescale" in keyword_parameters): 
        ximage=[int(i/lens[0])-1*(i<0) for i in x]
        yimage=[int(i/lens[1])-1*(i<0) for i in y]
        zimage=[int(i/lens[2])-1*(i<0) for i in z]
        x=[x[i]-ximage[i]*lens[0] for i in range(natom)]
        y=[y[i]-yimage[i]*lens[1] for i in range(natom)]
        z=[z[i]-zimage[i]*lens[2] for i in range(natom)]

    for i in range(len(x)):
#        print("i= ",i)
        outfile.write("%s %s %s %s %s %s %s "%(id[i],molecules[i],types[i],charges[i],x[i],y[i],z[i]))	
        if	("rescale" in keyword_parameters or "images" in keyword_parameters):
            outfile.write("%s %s %s\n"%(ximage[i],yimage[i],zimage[i]))
        else:
            outfile.write("\n") 
    outfile.write("\n")

    if vel_yes==1:
    
        outfile.write("Velocities\n\n")
        for i in range(len(velocities[0])):
            for j in range(len(velocities)):
                outfile.write("%s  "%velocities[j][i])
            outfile.write("\n")




    if nbonds>0:
        outfile.write("\nBonds\n\n")	
        for i in range(nbonds):
            outfile.write("%s %s %s %s\n"%(i+1,bonds[0][i],bonds[1][i],bonds[2][i]))

#	print angles

    if nangles>0:
        outfile.write("\nAngles\n\n")   
        for i in range(nangles):
            outfile.write("%s %s %s %s %s\n"%(i+1,angles[0][i],angles[1][i],angles[2][i],angles[3][i]))

    if ndihedrals>0:
        outfile.write("\nDihedrals\n\n")
        for i in range(ndihedrals):
            outfile.write("%s %s %s %s %s %s\n"%(i+1,dihedrals[0][i],dihedrals[1][i],dihedrals[2][i],dihedrals[3][i],dihedrals[4][i]))

    if nimpropers>0:
        outfile.write("\nImpropers\n\n")
        for i in range(nimpropers):
            outfile.write("%s %s %s %s %s %s\n"%(i+1,impropers[0][i],impropers[1][i],impropers[2][i],impropers[3][i],impropers[4][i]))

    outfile.close()
    return


