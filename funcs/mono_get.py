def	mono_get(file):
    import data_get
    infile=open(file,'r')
    mono_coor=data_get.data_get(file,'id',1)

    split_line=''
    while split_line!=['bv']:
        read_line=infile.readline()
        split_line=read_line.split()

    read_line=infile.readline()
    split_line=read_line.split()
    bvec=[]
    while split_line!=[]:
        bvec.append([float(i) for i in split_line])
        if len(bvec[-1])!=4:
            bvec[-1].append((sum([i**2 for i in bvec[-1]]))**.5)
        read_line=infile.readline()
        split_line=read_line.split()

    mono_coor.append(bvec)

    return mono_coor




def     mono_get_multi(file):
    import data_get

    print("")
    print("getting mono data")

    infile=open(file,'r')


    split_line=["sup"]

    type_found=0

    total_data=[]


    def	read_to(keys):
        key_found=0
        split_count=0
        while key_found==0:
            read_line=infile.readline()
            split_line=read_line.split()
            if len(split_line)>0:
                split_count=0
                for j in keys:
                    if split_line[0]==j:
                        key_found=1

            elif len(split_line)==0: split_count+=1
            if split_count>=100: key_found=1

        return split_line

#    print("post definition")
    split_line=read_to(["type","id"])
#    print( "split_line= ", split_line)
    type=int(split_line[1])
    id_count=0

    split_line=read_to(["type","id"])
#    print("split_line= ",split_line)
    while split_line!=[]:
        
        if split_line[0]=="type":
            type=int(split_line[1])
        elif split_line[0]=="id":
            id_count+=1
            mono_coor=data_get.data_get(file,'id',1)

            mono_coor=[[] for i in range(5)]
            sub_split=["hi"]
            read_line=infile.readline()
            sub_split=read_line.split()

            
            while sub_split!=[]:
                for i in range(5):
                    mono_coor[i].append(sub_split[i])
                read_line=infile.readline()
                sub_split=read_line.split()

            sub_split=''
            while sub_split!=['bv']:
                read_line=infile.readline()
                sub_split=read_line.split()

            read_line=infile.readline()
            sub_split=read_line.split()
            bvec=[]
            while sub_split!=[]:
                bvec.append([float(i) for i in sub_split])
                if len(bvec[-1])!=4:
                    bvec[-1].append((sum([i**2 for i in bvec[-1]]))**.5)
                read_line=infile.readline()
                sub_split=read_line.split()

            mono_coor.append(bvec)
            mono_coor.append(type)
            total_data.append(mono_coor)

        split_line=read_to(["type","id"])



    infile.close()
    infile=open(file,'r')
    split_line=read_to(["stereo"])
    read_line=infile.readline()
    split_line=read_line.split()
    relations=[[] for i in range(3)]
    while split_line !=[]:
        for i in range(3): relations[i].append(split_line[i])
        read_line=infile.readline()
        split_line=read_line.split()
    total_data.append(relations)

    return total_data	




