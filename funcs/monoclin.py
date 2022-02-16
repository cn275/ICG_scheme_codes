def monoclin(file,key ,*positional_parameters, **keyword_parameters):
    xyz=open(file)
    split_line=[]
    for i in range(2):
        read_line=xyz.readline()
        split_line=read_line.split()	
    time=[str(split_line[0])]
    for i in range(2):
        read_line=xyz.readline()
        split_line=read_line.split()
    natom=int(read_line)
    told=int(time[0])
    if time !=[str(key)]:	
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

    xy=0
    xz=0
    yz=0	
    xlen=(float(split_line[1])-float(split_line[0]))
    xy=float(split_line[2])

    read_line=xyz.readline()
    split_line=read_line.split()
    xz=float(split_line[2])
    ylen=(float(split_line[1])-float(split_line[0]))

    read_line=xyz.readline()
    split_line=read_line.split()
    yz=float(split_line[2])
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
        id=int(split_line[0])
        #		print id
        x[id-1]=float(split_line[2])*float(xlen)+xz*float(split_line[3])
        y[id-1]=float(split_line[3])*float(ylen)+xy*float(split_line[4])
        z[id-1]=float(split_line[4])*float(zlen)+yz*float(split_line[2])
        type[id-1]=int(split_line[1])
    xyz.close()
    return [x,y,z,type,[xlen,ylen,zlen,xy,xz,yz]]













