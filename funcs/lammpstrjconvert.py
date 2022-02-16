def lammpstrjconvert(x,y,z,id,type,xlen,ylen,zlen,filename,*positional_parameters,**keyword_parameters):
    import os
    
#    print("type= ", type)

    if ('append' in keyword_parameters): 
      if bool(keyword_parameters['append']):  outfile=open(filename,'a')
      else: outfile=open(filename,'w')

    else:
        outfile=open(filename,'w')

    if ('time' in keyword_parameters): t=int(keyword_parameters['time'])
    else:   t=0

    if ('dump' in keyword_parameters):
        if keywor_parameters['dump']=="True":
            dumpform=1



    outfile.write("ITEM: TIMESTEP\n")
    outfile.write("%s\n"%t)
    outfile.write("ITEM: NUMBER OF ATOMS\n")
    outfile.write(str(len(x)))
    outfile.write("\n")
    outfile.write("ITEM: BOX BOUNDS ") 



    if ("monoclin_tilts" in keyword_parameters): 
        print("soup")
        tilts=[float(i) for i in keyword_parameters["monoclin_tilts"]]
        outfile.write("xy xz yz ")
    else:
        tilts=[0.0 for i in range(3)]


    outfile.write(" pp pp pp")
    xy=tilts[0]
    xz=tilts[1]
    yz=tilts[2]
    outfile.write("\n")
    outfile.write("0 %s %s\n"%(xlen,xy))
    outfile.write("0 %s %s\n"%(ylen,xz))
    outfile.write("0 %s %s\n"%(zlen,yz))
    outfile.write("ITEM: ATOMS id type xs ys zs\n")

#    print len(id), len(type), len(x), len(y), len(z)

    if ('unscaled' in keyword_parameters):
        xlen=1
        ylen=1
        zlen=1

    for i in range(len(id)):
        outfile.write("%s %s %s %s %s\n"%(id[i],type[i],(x[i]-xz*y[i]/float(ylen))/float(xlen), (y[i]-xy*z[i]/float(zlen))/float(ylen), (z[i]-yz*x[i]/float(xlen))/float(zlen)))

    outfile.close()
    return
