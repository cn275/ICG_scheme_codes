def	mono_rot(refbv, inqbv,*positional_parameters,**keyword_parameters):
    import sys, math
    break_marker=0
    if len(refbv)!=len(inqbv):
        print("unequal numbers of bond vectors")
        break_marker=1
    if break_marker==1:
        sys.exit()
    # rescale the vectors if you would

    for i in range(len(refbv)):
        refbv[i]=[j/refbv[i][-1] for j in refbv[i]] 
    for i in range(len(inqbv)):
        inqbv[i]=[j/inqbv[i][-1] for j in inqbv[i]]

    print("refbv")
    for i in refbv: print(i)

    print("monobv")
    for i in inqbv: print(i)

    #we assume that refbv, and inqbv are to be matched element-to-element
    #Unless....
    if "element_match" in keyword_parameters:
        el_match=[int(i) for i in keyword_parameters["element_match"]]	
    else:
        el_match=[i for i in range(len(refbv))]

    if "mesh_density" in keyword_parameters:
        md=int(keyword_parameters["mesh_density"])
    else:
        md=10
    if "ncycles" in keyword_parameters:
        ncycle=int(keyword_parameters['ncycles'])
    else:
        ncycle=1


    angle_inc=3.1415*2.0/float(md)
    inq_rot=[[0 for i in range(3)] for j in range(3)]

    amin=0
    bmin=0
    dmin=0
    max_match=0
    max_id=0

    for cycle in range(ncycle):

        new_inc=angle_inc/(float(md)*.5)**cycle
        max_match=0
        max_id=0
        for i in range(md**3):
            a=new_inc*(i%md)+amin
            d=new_inc*(i/(md**2))+dmin
            b=new_inc*((i%(md**2))/md)+bmin

            ca=math.cos(a)
            cb=math.cos(b)
            cd=math.cos(d)
            sa=math.sin(a)
            sb=math.sin(b)
            sd=math.sin(d)


            print("angles= ",ca,cb,cd,sa,sb,sd)

            cur=0
            new_vec=[0,0,0]
            neg_dot=0
            for j in range(len(refbv)):
                bid=el_match[j]
                new_vec[0]=inqbv[bid][0]*ca*cb+inqbv[bid][1]*(ca*sb*sd-sa*cd)+inqbv[bid][2]*(ca*sb*cd+sa*sd)
                new_vec[1]=inqbv[bid][0]*sa*cb+inqbv[bid][1]*(sa*sb*sd+ca*cd)+inqbv[bid][2]*(sa*sb*cd-sa*sd*ca)
                new_vec[2]=-1*inqbv[bid][0]*sb+inqbv[bid][1]*(cb*sd)+inqbv[bid][2]*(cd*cb)
                dot=new_vec[0]*refbv[j][0]+new_vec[1]*refbv[j][1]+new_vec[2]*refbv[j][2]
                print("vec,dot= ",new_vec[0], new_vec[1], new_vec[2], dot)
                if dot<0: neg_dot=1
                cur+=dot
            print("a,b,d,dot= ",a,b,d,cur)
            if cur>max_match:
                min_angles=[a,b,d]
                max_match=cur
                max_id=i

		#redo the calculation but instead with a focused range
        amin=new_inc*(max_id%md-1)+amin
        bmin=new_inc*((max_id%(md**2))/md-1)+bmin	
        dmin=new_inc*(max_id/(md**2)-1)+dmin

        new_inc=angle_inc/float(md)*2.0

    if "dot_out" in keyword_parameters:
        output=max_match
    else:
        output=min_angles
    return	output
    







def	sup_mono_rot(refbv, inqbv,*positional_parameters,**keyword_parameters):

    import sys, math
    break_marker=0
    if len(refbv)!=len(inqbv):
        print("unequal numbers of bond vectors")
        break_marker=1
    if break_marker==1:
        sys.exit()





    # rescale the vectors if you would

    for i in range(len(refbv)):
        refbv[i]=[j/refbv[i][-1] for j in refbv[i]] 
    for i in range(len(inqbv)):
        inqbv[i]=[j/inqbv[i][-1] for j in inqbv[i]]


    #we assume that refbv, and inqbv are to be matched element-to-element
    #Unless....
    if "element_match" in keyword_parameters:
        el_match=[int(i) for i in keyword_parameters["element_match"]]	
    else:
        el_match=[i for i in range(len(refbv))]





    if "mesh_density" in keyword_parameters:
        md=int(keyword_parameters["mesh_density"])
    else:
        md=10
    if "ncycles" in keyword_parameters:
        ncycle=int(keyword_parameters['ncycles'])
    else:
        ncycle=1


    angle_inc=3.1415*2.0/float(md)
    inq_rot=[[0 for i in range(3)] for j in range(3)]


    amin=0
    bmin=0
    dmin=0
    max_match=0
    max_id=0









    for cycle in range(ncycle):

        new_inc=angle_inc/(float(md)*.5)**cycle
        #		print "angle_inc, new_inc= ",angle_inc, new_inc
        max_match=0
        max_id=0

        for i in range(md**3):
            a=new_inc*(i%md)+amin
            d=new_inc*(i/(md**2))+dmin
            b=new_inc*((i%(md**2))/md)+bmin

            #			print a,b,d	

            ca=math.cos(a)
            cb=math.cos(b)
            cd=math.cos(d)

            sa=math.sin(a)
            sb=math.sin(b)
            sd=math.sin(d)

            #			print "trig angles= ",ca,cb,cd,sa,sb,sd

            cur=0
            new_vec=[0,0,0]
            neg_dot=0
            for j in range(len(refbv)):
            #				print "ref bond vec= ",refbv[j]
			
                bid=el_match[j]
                new_vec[0]=inqbv[bid][0]*ca*cb+inqbv[bid][1]*(ca*sb*sd-sa*cd)+inqbv[bid][2]*(ca*sb*cd+sa*sd)
                new_vec[1]=inqbv[bid][0]*sa*cb+inqbv[bid][1]*(sa*sb*sd+ca*cd)+inqbv[bid][2]*(sa*sb*cd-sa*sd*ca)
                new_vec[2]=-1*inqbv[bid][0]*sb+inqbv[bid][1]*(cb*sd)+inqbv[bid][2]*(cd*cb)

                #				print "old vec, rot vec", inqbv[bid], new_vec	
                dot=new_vec[0]*refbv[j][0]+new_vec[1]*refbv[j][1]+new_vec[2]*refbv[j][2]
                #				print "dot= ", dot
                if dot<0: neg_dot=1
                # try enforcing that the bond vectors must be pointing in the same direction
                cur+=dot
            print("a,b,d,dot= ",a,b,d,cur)
            if cur>max_match:
                min_angles=[a,b,d]
                max_match=cur
                max_id=i

#		print "max_match= ",max_match
	

		#redo the calculation but instead with a focused range
        amin=new_inc*(max_id%md-1)+amin
        bmin=new_inc*((max_id%(md**2))/md-1)+bmin	
        dmin=new_inc*(max_id/(md**2)-1)+dmin

        new_inc=angle_inc/float(md)*2.0

#	print "max_match= ",max_match
#	print "min angles= ", min_angles

    if "dot_out" in keyword_parameters:
        output=max_match
    else:
        output=min_angles
    return	output
    







