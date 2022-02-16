def fileo(inputs,inputnames, filename,*positional_parameters,**keyword_parameters):
    import os
#	os.system("rm %s"%filename)
    if ('append' in keyword_parameters):
        if keyword_parameters['append']==True:
            output=open(filename,'a')
            output.write("\n\n")

    else:			
        output=open(filename,'w')


    if ("preamble" in keyword_parameters):
        for i in keyword_parameters["preamble"]:
            output.write("%s \n"%i)

#    output.write("\n")

    for j in range(len(inputnames)):
        output.write(inputnames[j])
        output.write(" ")
    output.write("\n")
    for i in range(len(inputs[0])):
        for j in range(len(inputs)):
            output.write(str(inputs[j][i]))
            output.write("	")
        output.write("\n")
    return 




def fileo_trans(inputs,inputnames, filename,*positional_parameters,**keyword_parameters):
    import os
#   os.system("rm %s"%filename)

    if ('append' in keyword_parameters):
        
        if keyword_parameters['append']==True:
            print("appending")
            output=open(filename,'a')
            output.write("\n\n")

    else:
        output=open(filename,'w')


    if ("preamble" in keyword_parameters):
        for i in keyword_parameters["preamble"]:
            output.write("%s \n"%i)

    output.write("\n\n")

    for j in range(len(inputnames)):
        output.write(inputnames[j])
        output.write(" ")
    output.write("\n")
    for i in range(len(inputs)):
        for j in range(len(inputs[i])):
            output.write(str(inputs[i][j]))
            output.write("  ")
        output.write("\n")

    output.close()
    return





def csv(inputs,inputnames, filename,*positional_parameters,**keyword_parameters):
    import os
#   os.system("rm %s"%filename)

    if ('append' in keyword_parameters):
        if keyword_parameters['append']=="True":
            output=open(filename,'a')
            output.write("\n\n")


    else:
        output=open(filename,'w')

    if ("header"  in keyword_parameters):
        if keyword_parameters["header"]==True:
            for j in range(len(inputnames)):
                output.write(inputnames[j])
                output.write(" ")
            output.write("\n")
    for i in range(len(inputs[0])):
        for j in range(len(inputs)):
            output.write(str(inputs[j][i]))
            if j+1!=len(inputs):
                output.write(",  ")
        output.write("\n")
    return
	
