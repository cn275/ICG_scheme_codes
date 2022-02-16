import math


def canberra_dist(v1,v2):
#    print("entering canberra_dist: ")
    sum_can=0

    for ni, i in enumerate(v1): 
#        print(i,v2[ni], math.fabs(i+v2[ni])/(math.fabs(i)+math.fabs(v2[ni])))
        if i==0 and v2[ni]==0:
#            print(i,v2[ni], 0)
            sum_can=sum_can+0
        else:
#            print(i,v2[ni], math.fabs(i-v2[ni])/(math.fabs(i)+math.fabs(v2[ni])))
            sum_can=sum_can+float(math.fabs(i-v2[ni]))/float(math.fabs(i)+math.fabs(v2[ni])) 
#    print("sum_can= ", sum_can)
    return  sum_can
