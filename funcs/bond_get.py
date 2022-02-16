def bond_get(file_name):		
# Here we read in the file and read and split the first line to initiate the loop

    file1 = open(str(file_name))
    read_line=file1.readline()
    split_line=read_line.split()
#    print('first split entry =',split_line[0])
    keynum=1
    # Here the Entry Index is the line at which we encounter the key phrase 
    # We assume that all Data files will some kind of keyword that symbolizes the
    # start of the actual data so we can ignore the rest
    entry_index=1
    keycount=0
    key="Bonds"
    # This block of commands is to help navigate through the file until we find the key
    # which is a keyword that signifies the next line is actual data
    # We read and split each line succesively looking for the first entry in a line to be key
    while keycount < keynum:
        while split_line == [] or split_line[0] != key:
            read_line=file1.readline()
            split_line=read_line.split()
            entry_index = entry_index + 1
        Nvar=len(split_line)
        if split_line[0] == key:
            keycount=keycount+1
            if keycount < keynum:
                read_line=file1.readline()
                split_line=read_line.split()
                entry_index = entry_index + 1
                read_line=file1.readline()
#                print(keycount)
#    print('entry_index=', entry_index)
    Nvar=4
    # These commands are now to create lists for each value
    # We simply keep appending the entries until we get to the end
    # The point of the len check is that if we encounter the phrase: "End Run", "Error"
    # or some other phrase which is shorter than the number of vars we
    # have then we can fill in those empty spots
    # with the 'empty'

    # the point of the phrase globals()['Var%s'%i] is so that we
    # can create the appropriate amount of lists to fix the issue
    # of a non static amount of variable we may encounter

    read_line=file1.readline()
    for i in range(Nvar):
        globals()['Var%s'%i]=[]

    while split_line != []:
        read_line=file1.readline()
        split_line=read_line.split()

        if split_line==[]:
#            print('read complete')
            v=2
        else:
            if len(split_line) < Nvar:
                for i in range(Nvar):
                    globals()['Var%s'%i].append('empty')

            else:
                for i in range(Nvar):
                    globals()['Var%s'%i].append(int(split_line[i]))
    output=[]
    for i in range(Nvar):
        output.append(globals()['Var%s'%i])                 
    return output

