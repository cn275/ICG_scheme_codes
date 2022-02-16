def mono_bond_get(file_name,key,keynum):
    print(file_name)
    print(key)
    print(keynum)
    keycount=0		
    # Here we read in the file and read and split the first line to initiate the loop
    file1 = open(str(file_name))
    read_line=file1.readline()
    split_line=read_line.split()

    while keycount < keynum:
        while split_line == [] or split_line[0] != key:
            read_line=file1.readline()
            split_line=read_line.split()



        Nvar=len(split_line)
        if split_line[0] == key:
            keycount=keycount+1
            if keycount < keynum:
                read_line=file1.readline()
                split_line=read_line.split()

	
    real_bonds=[[] for i in range(2)]
    connect_bonds=[[] for i in range(2)]
    while split_line != []:
        read_line=file1.readline()
        split_line=read_line.split()
        if split_line==[]:
            v=2
        else:
            if len(split_line[0].split('x'))==1 and len(split_line[1].split('x'))==1:
                real_bonds[0].append(split_line[0])
                real_bonds[1].append(split_line[1])
            else:
                connect_bonds[0].append(split_line[0])
                connect_bonds[1].append(split_line[1])




    return [real_bonds, connect_bonds]




def     mono_bond_get_multi(file):
    import data_get
    infile=open(file,'r')
    split_line=["sup"]
    type_found=0
    total_data=[]

    def     read_to(keys):
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








    split_line=read_to(["bonda"])
    print("split_line= ",split_line)
    while split_line!=[]:

        sub_split=["sup"]
        real_bonds=[[] for i in range(2)]
        connect_bonds=[[] for i in range(2)]

        read_line=infile.readline()
        sub_split=read_line.split()

        while sub_split != []:
            if len(sub_split[0].split('x'))==1 and len(sub_split[1].split('x'))==1:
                real_bonds[0].append(sub_split[0])
                real_bonds[1].append(sub_split[1])
            else:
                connect_bonds[0].append(sub_split[0])
                connect_bonds[1].append(sub_split[1])


#            print("bonds= ",[real_bonds,connect_bonds])

            read_line=infile.readline()
            sub_split=read_line.split()
        total_data.append([real_bonds,connect_bonds])
        print()
        print("total_data= ",total_data[-1])
        split_line=read_to(["bonda"])
        


    return total_data



		


