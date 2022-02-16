def data_get(file_name,key,keynum,*positional_parameters,**keyword_parameters):


    import num_check
    keycount=0		
    file1 = open(str(file_name))
    preamble=[]
    if ('preamble' in keyword_parameters):
        preadd=1	
    else:
        preadd=0

    if preadd==1:
        read_line=file1.readline()
        preamble.append(read_line)
    
    split_line=[]

    while keycount < keynum:
        read_line=file1.readline()
        split_line=read_line.split()
 #       print("split= ", split_line)
        if len(split_line)>0:
            if split_line[0]==key:
                keycount+=1

    if preadd==1: del preamble[-1]

#    print("vars= ", split_line)
    var_name=[i for i in split_line]
    #    print "var_name= ", var_name
    #    print len(var_name)


    variables=[[] for i in range(len(var_name))]
    Nvar=len(variables)
#    print("Nvar= ", Nvar)
    while split_line != []:
        read_line=file1.readline()
        split_line=read_line.split()
#        print("split_line= ", split_line)
        if split_line!=[]:
            skip=0
            if len(split_line) < Nvar:
                skip=1 
#                print("skip")
            else:
                for i in split_line:
                    if num_check.num_check(i) == 0:
                        skip=1

            if skip==0:
                for i in range(len(split_line)):
#                    print("i= ", i)
                    variables[i].append(split_line[i])

    output=[]
    for i in range(len(var_name)):
        output.append(variables[i])
    output.append(var_name)

    if preadd==1: output=[output,preamble]
    file1.close()
    return output



#------------------------------------------------------------------------------------------------------------------------------------------------------



def data_get_ncn(file_name,key,keynum,*positional_parameters,**keyword_parameters):


    import num_check
    keycount=0
    file1 = open(str(file_name))
    read_line=file1.readline()
    split_line=read_line.split()
    Nvar=len(split_line)
    for i in range(Nvar):
        globals()['Var%s'%i]=[]
        globals()['Var%s'%i].append(split_line[i])

    while split_line != []:
        read_line=file1.readline()
        split_line=read_line.split()

        if split_line==[]:
            v=2
#                        print 'read complete'
        else:
            skip=0
            if len(split_line) < Nvar:
                skip=1
            else:
                for i in split_line:
                    if num_check.num_check(i) == 0:
                        skip=1
            skip=0
            if skip==0:
                for i in range(Nvar):
                    globals()['Var%s'%i].append(split_line[i])

#here since we cant put the return command in a loop we need to create a loop to create the
# 'output' list that is a list of the lists we extracted               

    output=[]
    for i in range(Nvar):
        output.append(globals()['Var%s'%i])
    output.append([i+1 for i in range(Nvar)])

    file1.close()

    return output










def data_get_T(file_name,key,keynum):
    print(file_name)
    print(key)
    print(keynum)
    keycount=0
# Here we read in the file and read and split the first line to initiate the loop
    file1 = open(str(file_name))
    read_line=file1.readline()
    split_line=read_line.split()

#        print 'first split entry =',split_line[0]

# Here the Entry Index is the line at which we encounter the key phrase
# We assume that all Data files will some kind of keyword that symbolizes the
# start of the actual data so we can ignore the rest
    entry_index=1

# This block of commands is to help navigate through the file until we find the key
# which is a keyword that signifies the next line is actual data
# We read and split each line succesively looking for the first entry in a line to be key
    while keycount < keynum:
        while split_line == [] or split_line[0] != key:
            read_line=file1.readline()
            split_line=read_line.split()
            entry_index = entry_index + 1
#               print 'entry_index=', entry_index
        Nvar=len(split_line)
        if split_line[0] == key:
            keycount=keycount+1
            if keycount < keynum:
                read_line=file1.readline()
                split_line=read_line.split()
                entry_index = entry_index + 1
#       print keycount
#       print 'entry_index=', entry_index

    output=[]

    read_line=file1.readline()
    split_line=read_line.split()

    names=[]
    while split_line != []:

        output.append([split_line[i+1] for i in range(len(split_line)-1)])
        names.append(split_line[0])
        read_line=file1.readline()
        split_line=read_line.split()

    output.append(names)
#here since we cant put the return command in a loop we need to create a loop to create the
# 'output' list that is a list of the lists we extracted

    file1.close()
    return output












def data_get_single_line(file_name,key,keynum):
    print(file_name)
    print(key)
    print(keynum)
    keycount=0
    # Here we read in the file and read and split the first line to initiate the loop
    file1 = open(str(file_name))
    read_line=file1.readline()
    split_line=read_line.split()
    #        print 'first split entry =',split_line[0]

    # Here the Entry Index is the line at which we encounter the key phrase
    # We assume that all Data files will some kind of keyword that symbolizes the
    # start of the actual data so we can ignore the rest
    entry_index=1

# This block of commands is to help navigate through the file until we find the key
# which is a keyword that signifies the next line is actual data
# We read and split each line succesively looking for the first entry in a line to be key
    while keycount < keynum:
        while split_line == [] or split_line[0] != key:
            read_line=file1.readline()
            split_line=read_line.split()
            entry_index = entry_index + 1
#               print 'entry_index=', entry_index
        Nvar=len(split_line)
        if split_line[0] == key:
            keycount=keycount+1
            if keycount < keynum:
                read_line=file1.readline()
                split_line=read_line.split()
                entry_index = entry_index + 1
                #       print keycount
                #       print 'entry_index=', entry_index

#here since we cant put the return command in a loop we need to create a loop to create the
# 'output' list that is a list of the lists we extracted

    file1.close()
    return split_line








def data_get_string(file_name,key,keynum,*positional_parameters,**keyword_parameters):


    import num_check
    keycount=0
    # Here we read in the file and read and split the first line to initiate the loop
    file1 = open(str(file_name))

    preamble=[]
    if ('preamble' in keyword_parameters):
    #               print "getting preamble"
        preadd=1
    else:
        preadd=0

    read_line=file1.readline()
    if preadd==1:
        preamble.append(read_line)
    split_line=read_line.split()
    #       print split_line
    #        print 'first split entry =',split_line[0]

    # Here the Entry Index is the line at which we encounter the key phrase 
    # We assume that all Data files will some kind of keyword that symbolizes the
    # start of the actual data so we can ignore the rest
    entry_index=1

    # This block of commands is to help navigate through the file until we find the key
    # which is a keyword that signifies the next line is actual data
    # We read and split each line succesively looking for the first entry in a line to be key
    while keycount < keynum:
        while split_line == [] or split_line[0] != key:
            read_line=file1.readline()
            if preadd==1: preamble.append(read_line)
            split_line=read_line.split()
#            print "split_line= ", split_line
            entry_index = entry_index + 1
#              print 'entry_index=', entry_index      
        Nvar=len(split_line)
        if split_line[0] == key:
            keycount=keycount+1
            if keycount < keynum:
                read_line=file1.readline()
                if preadd==1: preamble.append(read_line)
                split_line=read_line.split()
#                print "split_line= ", split_line
                entry_index = entry_index + 1

    if preadd==1: del preamble[-1]
    #       print keycount
    #       print 'entry_index=', entry_index
    Nvar=len(split_line)
    var_name=split_line
    #       print "Nvar= ", Nvar
# These commands are now to create lists for each value
# We simply keep appending the entries until we get to the end
# The point of the len check is that if we encounter the phrase: "End Run", "Error"
# or some other phrase which is shorter than the number of vars we
# have then we can fill in those empty spots
# with the 'empty'

# the point of the phrase globals()['Var%s'%i] is so that we
# can create the appropriate amount of lists to fix the issue
# of a non static amount of variable we may encounter
    for i in range(Nvar):
        globals()['Var%s'%i]=[]

    while split_line != []:
        read_line=file1.readline()
        split_line=read_line.split()

        if split_line==[]:
            v=2
#                        print 'read complete'
        else:
            skip=0

            if len(split_line) < Nvar:
                skip=1

#            try: 
#                float(split_line[0])

#            except ValueError:
#                skip=1



            if skip==0:
                for i in range(Nvar):
                    globals()['Var%s'%i].append(split_line[i])

#here since we cant put the return command in a loop we need to create a loop to create the
# 'output' list that is a list of the lists we extracted               

    output=[]
    for i in range(Nvar):
        del globals()['Var%s'%i][-1]
        output.append(globals()['Var%s'%i])
    output.append(var_name)

    if preadd==1: output=[output,preamble]
    file1.close()
    return output
