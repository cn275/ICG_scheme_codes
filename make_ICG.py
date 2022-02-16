"""
All the import functions to execute the CG generation
'"""

ICGLIB="path to ICG lib of funcs. To be replaced with modules"
from collections import Counter
from functools import reduce
import  operator, math, itertools, string
import  os,sys, numpy as np
sys.path.append(ICGLIB)
import  pbcorrectsingle, avg_std, data_get, mono_get, coor_get, coor_get_lammps, canberra_dist
import  lammpstrjconvert, coor_get, mono_rot, gen_rot, bond_get, monoclin
import 	mono_bond_get, dotprod
import fileo
import  bond_get, angle_get, dihedral_get, lammpsinput
from inspect import currentframe, getframeinfo

import argparse
import time
import random
random.seed()



parser= argparse.ArgumentParser(description="Code for generating a cg scheme")
parser.add_argument("file",help="The lammps file which includes pertinent information")
parser.add_argument("cgfile",help="The lammps file which includes pertinent information")
parser.add_argument("-rings", help="Find the rings in a molecules to inform the coarse-graining scheme",action="store_true")

# max_schemes = limit the maximum number of output schemes
max_schemes=2000

args=parser.parse_args()

# convert masses of atoms into atom type
# rdkit may have something for me here
mass_dict={1:"H", 12:"C", 14:"N", 16:"O", 15:"C", 32:"S", 19:"F"}


#type number to generic letter
letter_dict={1:"A",2:"B",3:"C",4:"D",5:"E",6:"G"}

# easier way to get letter to from number
comp_codes=list(string.ascii_uppercase)

print("comp_codes= ", comp_codes)




#
#
#   SMALL FUNCTIONS
#
#
#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


def prod(iterable):
    return reduce(operator.mul, iterable, 1)



#assign each atom in the molecule to its group within the scheme
def get_group_id(n_aa,scheme):
    group_id=[-1 for i in range(n_aa)]
    for ni,i in enumerate(scheme):
        for j in i:
            group_id[j-1]=ni

    return group_id


#for each cg bead find what other cg beads it is bonded to
def     get_cg_bond_partners(scheme,aa_bond_partners, molsize):
    scheme_bond_partners=[[] for i in range(len(scheme))]
    scheme_id=[0 for i in range(molsize)]
    for ni, i in enumerate(scheme):
        for j in i:
            scheme_id[j-1]=ni+1
    for ni in range(molsize):
        sa=scheme_id[ni]
        for nk, k in enumerate(aa_bond_partners[ni]):
            sb=scheme_id[k-1]
            if sa>0 and sb>0 and sa!=sb:
                if scheme_bond_partners[sa-1].count(sb)==0:
                    scheme_bond_partners[sa-1].append(sb)
                if scheme_bond_partners[sb-1].count(sa)==0:
                    scheme_bond_partners[sb-1].append(sa)
    return scheme_bond_partners





#code to get the branches from the smile
def get_branches(smile):
    branches=[]
    pos=0
    while pos<len(smile):
        in_branch=1
        if smile[pos]=="(" and in_branch==0:
            branches.append(["("])
            branch_closed=1
            while branch_closed!=0:
                pos+=1
                branches[-1].append(smile[pos])
                if smile[pos]=="(": 
                    branch_closed=0
    return branches




def color_beads(n_aa,groups ):
	"""
	type the atoms by their group
	"""

	types=[0 for i in range(n_aa)]
	for n, i in enumerate(groups):
		for j in i:
			types[j-1]=n+1
	return types




def closed_str(check_str):
    closed= True
    stored=[]
    ii=0   
    while ii<len(check_str):
        si=check_str[ii]
        if si=="(":
            stored.append(si)
        elif si ==")":
            if len(stored)==0:
                closed=False
                break
            else:
                stored.pop()
        ii+=1 
    return closed and len(stored)==0



def contiguous_str(check_str):
    contiguous=True
    if closed_str(check_str):
        if check_str[0]=="(":
            #now we check if this is a singular branch
            contiguous=False
            needed_len=0
            for ii in range(len(check_str)):    
                if closed_str([j for nj, j in enumerate(check_str) if nj<=ii]): 
                    needed_len=ii+1
                    break
            if needed_len!=len(check_str):                     
                contiguous=False
    else:
        contiguous=False
    return  contiguous









# test closed str func
if 1==0:
    samp_str="(B(C))(C)"
    print("samp_str= ", samp_str)
    print("balanced?: ", closed_str(samp_str))
    print("contiquous?: ", contiguous_str(samp_str))


    sys.exit()



def get_mapped(n_aa,groups,cg_groups):
	"""
	Taking in the current CG groups and see which are mapped
	"""
	
	mapped=[0 for i in range(n_aa)]
	for n, i in enumerate(groups):
		for j in i:
			mapped[cg_groups.index(j)]=1
	return mapped


def get_next_possibles(n_aa,groups,aa_bond_partners,mapped,H_list,lead_group):

	"""
	Find which atoms will be added during the next step of the algorithm
	"""

	next_possibles=[0 for i in range(n_aa)]
	for j in groups[lead_group]:
		for k in aa_bond_partners[j-1]:
			if mapped[k-1]==0 and H_list[k-1]==1:
				next_possibles[k-1]=1                         
	return next_possibles


def get_terminal(n_aa,cg_group,aa_bond_partners,mapped,H_list):


    terminals=[ int(sum([H_list[i-1]  for i in  atomistic_bond_partners[j-1] if i in cg_group])>1)   for j in cg_group]
    for j in range(n_aa): 
        if mapped[j]==1: terminals[j]=1
    for nj, j in enumerate(mapped):
        if j==1:
            for nk, k in enumerate(aa_bond_partners[cg_group[nj]-1]):
                if k in cg_group:
                    terminals[cg_group.index(k)]=0        

    # Trying an additional piece of the loop
    for k in cg_group:
        if terminals[cg_group.index(k)]==1:
            for nj, j in enumerate(aa_bond_partners[k-1]):
                if H_list[j-1]==1 and (cg_group.count(j)==0):
                    terminals[cg_group.index(k)]=0



    return terminals



def make_colored_trj(n_aa, aa_coors, H_list, groups,lens,suffix):
	"""
	quickly make a colored file. Useful for troubleshooting
	"""
	types=color_beads(n_aa,groups)
	types=[j for i,j in enumerate(types) if H_list[i]==1]
	modded_coors=[[xyz for j, xyz in enumerate(aa_coors[i]) if H_list[j]==1] for i in range(3)] 
	lammpstrjconvert.lammpstrjconvert(modded_coors[0],modded_coors[1],modded_coors[2],[q+1 for q in range(len(modded_coors[0]))],types,lens[0],lens[1],lens[2],"test_%s.lammpstrj"%suffix)
	return


def get_group_done(n_aa,cg_group,aa_bond_partners,groups,H_list,dynamic_avg):
	"""
	Take in a mapping and see which are 'done'
	"""
	group_done=[int(len(i)>=dynamic_avg) for i in groups]
	print("groups= ", groups)
	print("init_group_done= ", group_done)
	print("cg_group= ", cg_group)
	mapped=get_mapped(n_aa,groups,cg_group)
	for ni, i in enumerate(groups):
		print("groups= ", i)
		if group_done[ni]==0:
			avail_adds=0
			for nj, j in enumerate(i):
				print("j= ", j)
				for nk, k in enumerate(aa_bond_partners[j-1]):
					if cg_group.count(k)!=0:
						print("k= ", k, cg_group.index(k))
						if mapped[cg_group.index(k)]==0 and H_list[k-1]==1:
							avail_adds+=1
			if avail_adds==0:
				group_done[ni]=1

	return  group_done



def start_find(n_aa,cg_group,aa_bond_partners,groups,H_list,dynamic_avg):
    mapped=get_mapped(n_aa,groups,cg_group)
    terminal=get_terminal(n_aa,cg_group,aa_bond_partners,mapped,H_list)
    group_done=get_group_done(n_aa,cg_group,aa_bond_partners,groups,H_list,dynamic_avg)
    return  mapped, terminal, group_done




"""
END OF FUNCTION DEFINITIONS


"""

sys.exit()


#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------



final_schemes=[]
#atomistic lfile
afile=args.file
#lfile for cg
cgfile=args.cgfile

# Read in atomistic coordinate
atomistic_coors=coor_get_lammps.coor_get_lammps_mol(afile,q="True")

#Read in atomistic bonds
atomistic_bonds=bond_get.bond_get(afile)

#Read in mass of each atom
atom_mass=coor_get_lammps.coor_get_lammps(afile,masses="True")[0]

#Get number of different atom types
ntypes=len(atom_mass)

#get the number of atoms
natom=len(atomistic_coors[0])

#get the molecules each atom belongs to
mols=[int(i) for i in atomistic_coors[5]]

#get the number of molecules
nmols=max(mols)

smile_backbone=[0 for i in range(len(atomistic_coors[0]))]
copy_smile_backbone=[0 for i in range(len(atomistic_coors[0]))]

# get molecule size
molsize=int(len(atomistic_coors[0])/nmols)

# get number of bonds in a molecule
bond_per_mol=len(atomistic_bonds[0])/nmols

#reduce the coodinates to that of a single molecules
for i in range(len(atomistic_coors)-1): atomistic_coors[i]=[j for nj, j in enumerate(atomistic_coors[i]) if nj<molsize]

#same for bonds
for i in range(3):  atomistic_bonds[i]=[j for nj, j in enumerate(atomistic_bonds[i]) if nj<bond_per_mol]

#same for cg as we did for aa above
cg_coors=coor_get_lammps.coor_get_lammps_mol(cgfile,q="True")
cg_bonds=bond_get.bond_get(cgfile)
n_cg=len(cg_coors[0])/nmols
bond_per_cg=len(cg_bonds[0])/nmols
n_cg_types=max([int(i) for i in cg_coors[3]])


#get a count of each type in the types
# new way to do it:
# a=Counter(type)
cg_type_count=[len([i for i in cg_coors[3] if int(i)==j+1]) for j in range(n_cg_types)]
aa_type_count=[len([i for i in atomistic_coors[3] if int(i)==j+1]) for j in range(ntypes)]

#reduce cg coors to single molecule
for i in range(len(cg_coors)-1): cg_coors[i]=[j for nj, j in enumerate(cg_coors[i]) if nj<n_cg]
for i in range(3):  cg_bonds[i]=[j for nj, j in enumerate(cg_bonds[i]) if nj<bond_per_cg]
cg_molsize=len(cg_coors[0])


# get the separation distance map
# can redefine using numpy array initialization
cg_sep_map=[[[cg_molsize+1 for i in range(cg_molsize)] for j in range(cg_molsize)] for k in range(nmols)]

#By definition the self sep distance is 0
for i in range(nmols):
    for j in range(cg_molsize):
        cg_sep_map[i][j][j]=0


#initialize cg bond partners
cg_bond_partners=[[] for i in range(cg_molsize)]

# calculate the bond partners
for i in range(len(cg_bonds[2])):
    ba=int(cg_bonds[2][i])
    bb=int(cg_bonds[3][i])
    bond_mol=1
    baind=ba-1
    bbind=bb-1
    cg_sep_map[0][baind][bbind]=1
    cg_sep_map[0][bbind][baind]=1
    cg_bond_partners[baind].append(bbind+1)
    cg_bond_partners[bbind].append(baind+1)


# calculate the separation map
for i in range(nmols):
    for j in range(cg_molsize):
        level=j+1
        for k in range(cg_molsize):
            for l in range(cg_molsize):
                if cg_sep_map[i][k][l]==level:
                    for m in range(cg_molsize):
                        if cg_sep_map[i][l][m]==1 and m>k:
                            if cg_sep_map[i][k][m]>level+1:
                                cg_sep_map[i][k][m]=level+1
                                cg_sep_map[i][m][k]=level+1


#find where there is branching
cg_branch_point=[int(len([i for i in cg_sep_map[0][j] if i==1])>2) for j in range(cg_molsize)]

#find your branch partners
cg_branch_partners=[[ni+1 for ni, i in enumerate(cg_sep_map[0][j]) if i==1] for j in range(cg_molsize)]

# Determine which atoms are H or not
type_dict=["" for i in atom_mass]
for i in range(len(type_dict)):
	if int(round(atom_mass[i],0))==1:
		type_dict[i]="H"
	else:
		type_dict[i]="NH"



#create a dictionary of atom type and wether it is H or not
type_dict=dict((i+1,type_dict[i]) for i in range(len(type_dict)))

# get atoms that aren't H
H_list=[1-int(type_dict[int(i)]=="H") for i in atomistic_coors[3]]

#get sim box dimensions
lens=[i for i in atomistic_coors[-1]]

#get number of atoms and number of nonH atoms
n_aa=len(atomistic_coors[0])
n_nah=sum(H_list)


#new atomistic coordinates
new_atomistic=[[0 for i in range(len(atomistic_coors[0]))] for j in range(3)]

#get a bonds
abonds=bond_get.bond_get(afile)

#get atom bond partners
atomistic_bond_partners=[[] for i in range(len(atomistic_coors[0]))]
for i in range(len(atomistic_bonds[0])):
	atomistic_bond_partners[atomistic_bonds[2][i]-1].append(atomistic_bonds[3][i])
	atomistic_bond_partners[atomistic_bonds[3][i]-1].append(atomistic_bonds[2][i])

aa_bond_partners=[[j for j in i] for i in atomistic_bond_partners]

if nmols>1:
    print("Can't have more than one molecule right now")
    sys.exit()

#initialize the basic info on the atomistic
mol_size=[sum([1 for i in mols if i==(j+1)]) for j in range(nmols)]
mapped=[0 for i in range(n_aa)]
atom_in_mol=[[i+1 for i in range(natom) if mols[i]==j+1] for j in range(nmols)]
nh_atom_in_mol=[[i+1 for i in range(natom) if mols[i]==j+1 and H_list[i]==1] for j in range(nmols)]
id_in_mol=[0 for i in range(natom)]


#assign each atom an id within the molecule
for i in atom_in_mol:
        mol_id_count=0
        for j in i:
                id_in_mol[j-1]=mol_id_count
                mol_id_count+=1

#generate the molecular sep map
mol_sep_map=[[[molsize+1 for j in range(len(atom_in_mol[k]))] for i in range(len(atom_in_mol[k]))] for k in range(nmols)]
mol_natom=[len(atom_in_mol[i]) for i in range(nmols)]

#calculate the sep maps
for j in range(len(atomistic_bonds[0])):
    a1=int(atomistic_bonds[2][j])
    a2=int(atomistic_bonds[3][j])
    mol_of_bond=mols[a1-1]
    a1_id=atom_in_mol[mol_of_bond-1].index(a1)
    a2_id=atom_in_mol[mol_of_bond-1].index(a2)
    mol_sep_map[mol_of_bond-1][a1_id][a2_id]=1
    mol_sep_map[mol_of_bond-1][a2_id][a1_id]=1


for i in range(nmols):	
    for j in range(int(molsize)):
        level=j+1
        for k in range(mol_natom[i]):
            katm=atom_in_mol[i][k]
            for l in range(mol_natom[i]):
                latm=atom_in_mol[i][l]
                if mol_sep_map[i][k][l]==level:
                    for m in range(mol_natom[i]):
                        matm=atom_in_mol[i][m]
                        if mol_sep_map[i][l][m]==1 and m>k:
                            if level+1<mol_sep_map[i][k][m]:
                                mol_sep_map[i][k][m]=level+1
                                mol_sep_map[i][m][k]=level+1
    for j in range(mol_natom[i]):
        for k in range(mol_natom[i]):
            if j>k:
                mol_sep_map[i][j][k]=mol_sep_map[i][k][j]


#detemrine which nonH atoms are terminal
terminal=[ int(sum([H_list[i-1]  for i in  atomistic_bond_partners[j]])>1)   for j in range(n_aa)]

#output colored molecule showing the terminal atoms
lammpstrjconvert.lammpstrjconvert(atomistic_coors[0],atomistic_coors[1],atomistic_coors[2],[i+1 for i in range(n_aa)],terminal,lens[0],lens[1],lens[2],"terminal_test.lammpstrj")



for i in range(0):
	mol=mols[i]
	for j in range(natom):
		termination_map[j][i]=sum([ int(terminal[k-1]==0) for k in atom_in_mol[mol-1] if (mol_sep_map[mol-1][id_in_mol[j]][id_in_mol[k-1]]==i+1)])

all_groups=[]
temp=[]
n_mapped=0
orig_terminal=[terminal[i]*H_list[i] for i in range(n_aa)]


#get first degree neighs
d_i=[sum([ int(j==1) for nj, j in enumerate(mol_sep_map[0][i]) if H_list[nj]==1]) for i in range(natom)]

c_i=[0 for i in range(natom)]

dN_i=[sum([ int(j==2) for j in mol_sep_map[0][i]])/d_i[i] for i in range(natom)]
c_Ni=[0 for i in range(natom)]
egos=[[j+1 for j in range(natom) if mol_sep_map[0][i][j]==1 and H_list[j]==1] for i in range(natom)]
E_egoi=[sum([len([1 for k in range(natom) if mol_sep_map[0][j-1][k]==1 and egos[i].count(k+1)==1 and k!=i ]) for j in egos[i]])/2 for i in range(natom)]
E0_egoi=[sum([len([1 for k in range(natom) if mol_sep_map[0][j-1][k]==1 and egos[i].count(k+1)==0 and k!=i ]) for j in egos[i]]) for i in range(natom)]
ego_N_neigh=[0 for i in range(natom)]

for i in range(natom):
    ego_neigh=[]
    for j in egos[i]:
        for k in egos[j-1]:
            if ego_neigh.count(k)==0 and k-1!=i:
                ego_neigh.append(k)
    ego_N_neigh[i]=len(ego_neigh)

feature_vec=[[d_i[i],c_i[i],dN_i[i],c_Ni[i],E_egoi[i], E0_egoi[i], ego_N_neigh[i]] for i in range(natom)]



#--------------------------------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------------------------------



branch_point=[int(len([i for i in mol_sep_map[0][j] if i==1])>2) for j in range(molsize)]
branch_partners=[[ni+1 for ni, i in enumerate(mol_sep_map[0][j]) if i==1] for j in range(molsize)] 
smiles=[]
branch_depths=[]
cg_smiles=[]

smiles_id=[]
cg_smiles_id=[]

atomistic_coors[3]=[int(i) for i in atomistic_coors[3]]
cg_coors[3]=[int(i) for i in cg_coors[3]]
smile_type=[mass_dict[int(round(atom_mass[i-1],0))] for i in atomistic_coors[3]]

print("smile_type= ", smile_type)
print(H_list)

cg_smile_type=[comp_codes[i-1] for i in cg_coors[3]]
max_type=max([max(atomistic_coors[3]),max(atomistic_coors[3])])
for i in smile_type: 
    if comp_codes.count(i)>0:
        comp_codes.remove(i)



dynamic_avg=int(n_nah/len(cg_coors[0]))
print("dynamic_avg= ", dynamic_avg)


#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

def get_group_term_map(group1,aa_bond_partners, max_level):
    termination_map=[[0 for i in range((max_level))] for j in group1]
    sub_group_terminal=[int( sum([k for nk, k in enumerate(mol_sep_map[0][j-1]) if k==1])==1 ) for j in group1]

    extended_group=[i for i in group1]
    for i in group1:
        for j in range(molsize):
            if mol_sep_map[0][i-1][j]==1 and terminal[j]==1:
                if extended_group.count(j+1):
                    extended_group.append(j+1)

    for i in range(max_level):
        for nj, j in enumerate(group1):

            termination_map[nj][i]=sum([ terminal[k-1] for nk,k  in enumerate(extended_group) if (mol_sep_map[0][j-1][k-1]==i+1)])

    return  termination_map

#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


def compare_aa_groups(group1,group2, *positional_parameters, **keyword_parameters):


#    print("groups= ", group1, group2)

    swapped=0
    group1_copy=[i for i in group1]
    group2_copy=[i for i in group2]

    if len(group1)>len(group2):
        swapped=1
        group1=[i for i in group2_copy]
        group2=[i for i in group1_copy]

#    print("groups_next= ", group1, group2)        


    differing_dicts=[]
    n_cg=len(group1)
    n_cg2=len(group2)

    bond_out_1=[int(len([j for j in aa_bond_partners[i-1] if group1.count(j)==0 and H_list[j-1]==1])>0) for i in group1]
    bond_out_2=[int(len([j for j in aa_bond_partners[i-1] if group2.count(j)==0 and H_list[j-1]==1])>0) for i in group2]


    max_level1=max([max([mol_sep_map[0][i-1][j-1] for j in group1]) for i in group1])
    max_level2=max([max([mol_sep_map[0][i-1][j-1] for j in group2]) for i in group2])

    group1_term_map=get_group_term_map(group1, aa_bond_partners, max_level1)
    group2_term_map=get_group_term_map(group2, aa_bond_partners, max_level2)

    group1_sep_map=[[mol_sep_map[0][i-1][j-1] for j in group1] for i in group1]
    group2_sep_map=[[mol_sep_map[0][i-1][j-1] for j in group2] for i in group2]

    possibility_matrix=[[1 for i in range(n_cg2)] for j in range(n_cg)]


    if ("known_pairs" in keyword_parameters):
        for i in keyword_parameters["known_pairs"]:
            if swapped==0:
                i1=group1.index(i[0])
                i2=group2.index(i[1])
            else:
                i1=group1.index(i[1])
                i2=group2.index(i[0])
            possibility_matrix[i1]=[int(j==i2) for j in range(n_cg2)]
            for j in range(n_cg): possibility_matrix[j][i2]=int(j==i1)
            for j in range(n_cg):
                for k in range(n_cg2):
                        if mol_sep_map[0][i[0]-1][group1[j]-1]!=mol_sep_map[0][i[1]-1][group2[k]-1]:
                            possibility_matrix[j][k]=0

    def single_fit_elim():
        for j in range(n_cg):
            if sum(possibility_matrix[j])==1:
                one_fit=possibility_matrix[j].index(1)
                for k in range(n_cg):
                    for l in range(n_cg2):
                        if mol_sep_map[0][group1[j]-1][group1[k]-1]!=mol_sep_map[0][group2[one_fit]-1][group2[l]-1]:
                            possibility_matrix[k][l]=0

#        print("after elim possibility_matrix")
#        for i in possibility_matrix: print(i)
        return

    single_fit_elim()

    for j in range(n_cg):
        for k in range(n_cg2):
            if possibility_matrix[j][k]==1:
                for l in range(max_level1):
                    if group1_term_map[j][l]>group2_term_map[k][l]:
                        possibility_matrix[j][k]=0

    single_fit_elim()
    single_fit_elim()

#    print("after term map possibility_matrix")
#    for i in possibility_matrix: print(i)



    def select_degeneracies(possibility_matrix):
        deg1=0
        while deg1!=n_cg:
            deg_ids=[]
            deg1=n_cg
            for i in range(len(possibility_matrix)):
                if sum(possibility_matrix[i])>1:
                    deg1=i
                    deg_ids=[j for j in range(n_cg) if possibility_matrix[deg1][j]==1]
                    break

#            print("deg1= ", deg1)
#            print("deg_ids= ", deg_ids)

            for ni, i in enumerate(deg_ids):
                if ni==0:
                    for j in range(n_cg):
                        if j!=deg_ids[ni]:
                            possibility_matrix[deg1][j]=0
#                else:
#                    possibility_matrix[i][deg1]=0

#            print("After selecting for degen")                
            sum_pos=0
            for j in possibility_matrix: 
#                print(j)
                if sum(j)==0: sum_pos+=1
    
            if sum_pos>0: 
                print("we quitting")
                sys.exit()

            single_fit_elim()

#            print("After elim for degen") 
#            for j in possibility_matrix: 
#                print(j)

        return

#    print("before deg possibility_matrix")
#    for i in possibility_matrix: print(i)    
    select_degeneracies(possibility_matrix)
    cant_fit=0
#    print("final possibility_matrix")
#    for i in possibility_matrix: print(i)


    for i in possibility_matrix:
        if sum(i)==0:
            cant_fit=1


    done=1
    for i in possibility_matrix:
        if sum(i)!=1:
            done=0

    if done==1:
        if swapped==0:
            differing_dicts.append(dict((group1[i],group2[[possibility_matrix[j][i] for j in range(n_cg)].index(1)]) for i in range(n_cg)))
        else:
            differing_dicts.append(dict((group2[[possibility_matrix[i][j] for j in range(n_cg)].index(1)],group1[i]) for i in range(n_cg)))

    elif cant_fit==1:
        v=2

    else:
        print("sheeet not done")
        sys.exit()

#    print("differing_dict= ", differing_dicts)
    return differing_dicts










def     find_rings():
    in_ring=[0 for i in range(natom)]
    rings=[[] for i in range(nmols)]
    for q in range(nmols):
        group=[]
        q_mol_size=len(atom_in_mol[q])
        for i in atom_in_mol[q]:
            base=i
            if in_ring[i-1]==0 and terminal[i-1]!=0:
#                print()
#                print("i= ", i)
                ogb=[0 for j in range(q_mol_size)]
                lays=[0 for j in range(q_mol_size)]
                lays[id_in_mol[i-1]]=-1
                bnum=1
                for j in atomistic_bond_partners[i-1]:
                    if terminal[j-1]!=0:
                        ogb[id_in_mol[j-1]]=bnum
                        lays[id_in_mol[j-1]]=1
                        bnum+=1

                for j in range(2):
                    for k in range(q_mol_size):
                        if lays[k]==1+j and terminal[atom_in_mol[q][k]-1]!=0:
                            for l in atomistic_bond_partners[k]:
                                if terminal[l-1]!=0 and lays[id_in_mol[l-1]]==0:
                                    lays[id_in_mol[l-1]]=j+2
                                    ogb[id_in_mol[l-1]]=ogb[k]

#                print("lays= ", lays)
                

                
                for j in range(q_mol_size):
                    if lays[j]>0:
                        for k in atomistic_bond_partners[atom_in_mol[q][j]-1]:
                            if lays[id_in_mol[k-1]]>0:
                                if ogb[id_in_mol[k-1]]!=ogb[j]:
                                    loop_size=(lays[j]+lays[id_in_mol[k-1]])+1
#                                    print("loop_size= ", loop_size)
                                    if loop_size%2==1:
                                        group=[base, atom_in_mol[q][j], k]
                                        lay_end_1=lays[id_in_mol[group[1]-1]]
                                        lay_end_2=lays[id_in_mol[group[2]-1]]
#                                        print("group= ", group)
#                                        print("lay_end_1= ", lay_end_1)
#                                        print("lay_end_2= ", lay_end_2)
                                        start=group[1]
                                        while lay_end_1!=1:
                                            for m in atomistic_bond_partners[start-1]:
                                                if lays[id_in_mol[m-1]]==lay_end_1-1 and ogb[m-1]==ogb[start-1]:
                                                    group.append(m)
                                                    lay_end_1-=1
                                                    start=m
                                        start=group[2]
                                        while lay_end_2!=1:
                                            for m in atomistic_bond_partners[start-1]:
                                                if lays[id_in_mol[m-1]]==lay_end_2-1 and ogb[id_in_mol[m-1]]==ogb[id_in_mol[start-1]]:
                                                    group.append(m)
                                                    lay_end_2-=1
                                                    start=m

                                    if loop_size%2==0:
                                        group=[base, atom_in_mol[q][j], k]
                                        lay_end_1=lays[id_in_mol[group[1]-1]]
                                        lay_end_2=lays[id_in_mol[group[2]-1]]
                                        start=group[1]
                                        while lay_end_1!=1:
                                            for m in atomistic_bond_partners[start-1]:
                                                if lays[id_in_mol[m-1]]==lay_end_1-1 and ogb[m-1]==ogb[start-1]:
                                                    group.append(m)
                                                    lay_end_1-=1
                                                    start=m
                                        start=group[2]
                                        while lay_end_2!=1:
                                            for m in atomistic_bond_partners[start-1]:
                                                if lays[id_in_mol[m-1]]==lay_end_2-1 and ogb[id_in_mol[m-1]]==ogb[id_in_mol[start-1]]:
                                                    group.append(m)
                                                    lay_end_2-=1
                                                    start=m
                    if len(group)>0:
                        fail=0
                        for v in rings[q]:
                            if fail==0:
                                for w in group:
                                    if v.count(w)==1:
                                        fail=1
                        if fail==0:
                            rings[q].append(group)
                            for w in group: in_ring[w-1]=1
                            break
    return rings



def get_cg_smile(atom_to_use, start_str, *positional_parameters, **keyword_parameters):
    tempsmile=''
    first_end=0
    if ("pin_end" in keyword_parameters): first_end=int(keyword_parameters["pin_end"])
    if len(atom_to_use)>1:
        if first_end==0:
            largest_sep=0
            sep_pair=[]
            for i in range(cg_molsize):
                for j in range(cg_molsize):
                    if j>i:
                        cur_sep= cg_sep_map[0][i][j]
                        if cur_sep==largest_sep: sep_pair.append([i,j])
                        if cur_sep>=largest_sep:
                            sep_pair=[[i,j]]
                            largest_sep=cur_sep
        else:
            largest_sep=0
            sep_pair=[]
            for j in range(molsize):
                if atom_to_use.count(j+1)==1 and j!=first_end-1:
                    cur_sep= cg_sep_map[0][first_end-1][j]
                    if cur_sep==largest_sep:
                        sep_pair.append([first_end-1,j])
                    if cur_sep>=largest_sep:
                        sep_pair=[[first_end-1,j]]
                        largest_sep=cur_sep
        if first_end==0: get_cg_smile(atom_to_use,'',pin_end=sep_pair[0][1]+1)
        sep_path=[0 for i in range(largest_sep+1)]
        sep_path[0]=sep_pair[0][0]
        sep_path[-1]=sep_pair[0][-1]
        for i in range(largest_sep-1):
            for j in cg_bond_partners[sep_path[i]]:
                if atom_to_use.count(j)==1:
                    if cg_sep_map[0][j-1][sep_path[-1]]==largest_sep-1-i: 
                        sep_path[i+1]=j-1
        tempsmile+=cg_smile_type[sep_pair[0][0]]
        for j in cg_branch_partners[sep_path[0]]:
            if sep_path.count(j-1)==0 and atom_to_use.count(j)==1:
                send_group=[j]+[k+1 for k in range(cg_molsize) if cg_sep_map[0][k][j-1]<cg_sep_map[0][k][sep_path[i+1]] and k!=sep_path[i+1] and k+1!=j]
                tempsmile+="(%s)"%get_cg_smile(send_group,'',pin_end=j)
        for i in range(largest_sep):
            if cg_branch_point[sep_path[i+1]]==0:
                tempsmile+=cg_smile_type[sep_path[i+1]]
            else:
                tempsmile+=cg_smile_type[sep_path[i+1]]
                for j in cg_branch_partners[sep_path[i+1]]:
                    if sep_path.count(j-1)==0 and atom_to_use.count(j)==1:
                        send_group=[j]+[k+1 for k in range(cg_molsize) if cg_sep_map[0][k][j-1]<cg_sep_map[0][k][sep_path[i+1]] and k!=sep_path[i+1] and k+1!=j]
                        tempsmile+="(%s)"%get_cg_smile(send_group,'',pin_end=j)
    else:
        for i in atom_to_use: tempsmile+=cg_smile_type[i-1]
    if len(atom_to_use)==cg_molsize:
        cg_smiles.append([k for k in tempsmile])
    return tempsmile


get_cg_smile([i+1 for i in range(cg_molsize)],'')







overall_smile_backbone=[]

def get_smile(atom_to_use, start_str, depth,*positional_parameters, **keyword_parameters):
#    print('atom to use= ', atom_to_use)

    smile_backbone=[0 for i in range(molsize)]

    tempsmile=''
    temp_branches=[]
    temp_id=[]
    first_end=0
    second_end=0
    if ("pin_end" in keyword_parameters):
        first_end=int(keyword_parameters["pin_end"])
    if ("pin_second" in keyword_parameters):
        second_end=int(keyword_parameters["pin_second"])

    if len(atom_to_use)>1:
        if first_end==0:
            #find largest sep
            largest_sep=0
            sep_pair=[]
            for i in range(molsize):
                if H_list[i]==1:
                    for j in range(molsize):
                        if H_list[j]==1:
                            if j>i:
                                cur_sep= mol_sep_map[0][i][j]
                                if cur_sep==largest_sep:
                                    sep_pair.append([i,j]) 
                                if cur_sep>largest_sep: 
                                    sep_pair=[[i,j]]
                                    largest_sep=cur_sep

                                 
        else:
            largest_sep=0
            sep_pair=[]
            for j in range(molsize):
                if H_list[j]==1:
                    if atom_to_use.count(j+1)==1 and j!=first_end-1:
                        cur_sep= mol_sep_map[0][first_end-1][j]
                        if cur_sep==largest_sep:
                            sep_pair.append([first_end-1,j]) 
                        if cur_sep>=largest_sep: 
                            sep_pair=[[first_end-1,j]]
                            largest_sep=cur_sep

        if first_end==0:
            get_smile(atom_to_use,'',depth,pin_end=sep_pair[0][1]+1, pin_second=sep_pair[0][0]+1)

        for nj, j in enumerate(sep_pair):
            if depth==0:
                smile_backbone[sep_pair[nj][1]]=1
                smile_backbone[sep_pair[nj][0]]=1
#                print("in backbone 1= ", sep_pair[nj][1], sep_pair[nj][0])
            if nj!=0:
                get_smile(atom_to_use,'',depth,pin_end=sep_pair[nj][1]+1, pin_second=sep_pair[nj][0]+1)
                get_smile(atom_to_use,'',depth,pin_end=sep_pair[nj][0]+1, pin_second=sep_pair[nj][1]+1)


        sep_path=[0 for i in range(largest_sep+1)]
        sep_path[0]=sep_pair[0][0]


        if second_end!=0:
            sep_path[-1]=second_end-1
        else:
            sep_path[-1]=sep_pair[0][-1]


        for i in range(largest_sep-1):
            for j in atomistic_bond_partners[sep_path[i]]:
                if atom_to_use.count(j)==1 and H_list[j-1]==1:
                    if mol_sep_map[0][j-1][sep_path[-1]]==largest_sep-1-i:
                        sep_path[i+1]=j-1

        tempsmile+=smile_type[sep_pair[0][0]]
        temp_id.append(sep_pair[0][0]+1)

        for j in branch_partners[sep_path[0]]:
            if sep_path.count(j-1)==0 and atom_to_use.count(j)==1 and H_list[j-1]==1:
                send_group=[j]+[k+1 for k in range(molsize) if mol_sep_map[0][k][j-1]<mol_sep_map[0][k][sep_path[i+1]] and k!=sep_path[i+1] and k+1!=j and H_list[k]==1]
                print("send_group= ", send_group)
                sub_smile=get_smile(send_group,'',depth+1,pin_end=j)
                tempsmile+="(%s)"%sub_smile[0]
                temp_id.append("(")
                temp_id.extend([k for k in sub_smile[1]])
                temp_id.append(")")

        for i in range(largest_sep):
            if depth==0:
#                print("in backbone 3= ", sep_path[i+1])
                smile_backbone[sep_path[i+1]]=1
            if branch_point[sep_path[i+1]]==0:
                tempsmile+=smile_type[sep_path[i+1]]     
                temp_id.append(sep_path[i+1]+1)

            else:
                tempsmile+=smile_type[sep_path[i+1]]
                temp_id.append(sep_path[i+1]+1)
                for j in branch_partners[sep_path[i+1]]:
                    if sep_path.count(j-1)==0 and atom_to_use.count(j)==1 and  H_list[j-1]==1:
                        send_group=[j]+[k+1 for k in range(molsize) if mol_sep_map[0][k][j-1]<mol_sep_map[0][k][sep_path[i+1]] and k!=sep_path[i+1] and k+1!=j and H_list[k]==1]
                        print("send_group= ", send_group)
                        sub_smile=get_smile(send_group,'',depth+1,pin_end=j)
                        tempsmile+="(%s)"%sub_smile[0]
                        temp_id.append("(")
                        temp_id.extend([k for k in sub_smile[1]])
                        temp_id.append(")") 
 
    else:
        for i in atom_to_use: 
            tempsmile+=smile_type[i-1] 
            temp_id.append(i)

    if depth==0:
        smiles.append([k for k in tempsmile])
        smiles_id.append([k for k in temp_id])
        overall_smile_backbone.append(smile_backbone)

    return [tempsmile, temp_id, depth]













copy_overall_smile_backbone=[]




def get_smile_rings(atom_to_use, start_str, depth, n_rings, keep_in_mol, *positional_parameters, **keyword_parameters):


#    print("atoms_to_use= ", atom_to_use)


    copy_smile_backbone=[0 for i in range(molsize+n_rings)]

    tempsmile=''
    temp_branches=[]
    temp_id=[]
    first_end=-1
    second_end=0
    if ("pin_end" in keyword_parameters):
        first_end=int(keyword_parameters["pin_end"])
#        print("first_end= ", first_end)
    if ("pin_second" in keyword_parameters):
        second_end=int(keyword_parameters["pin_second"])

    if len(atom_to_use)>1:
        if first_end==-1:
            #find largest sep
            largest_sep=0
            sep_pair=[]
            for i in range(molsize+n_rings):
                if atom_to_use.count(i+1)==1 and copy_H_list[i]==1:
#                    print("checking for sep pair of: ", i+1)
                    for j in range(molsize+n_rings):
                        if atom_to_use.count(j+1)==1 and j>i and copy_H_list[j]==1:
                            cur_sep= copy_mol_sep_map[0][i][j]
                            if cur_sep==largest_sep:
                                sep_pair.append([i,j])
                            if cur_sep>largest_sep:
                                sep_pair=[[i,j]]
                                largest_sep=cur_sep

        else:
            largest_sep=0
            sep_pair=[]
            for j in range(molsize+n_rings):
#                print("j= ",j, atom_to_use.count(j+1))
                if atom_to_use.count(j+1)==1 and j!=first_end-1 and copy_H_list[j]==1:
                    cur_sep= copy_mol_sep_map[0][first_end-1][j]
                    if cur_sep==largest_sep:
                        sep_pair.append([first_end-1,j])
                    if cur_sep>=largest_sep:
                        sep_pair=[[first_end-1,j]]
                        largest_sep=cur_sep

#        print("sep_pair= ", sep_pair[0], depth)

        if first_end==-1:
            get_smile_rings(atom_to_use,'',depth,n_rings, keep_in_mol,pin_end=sep_pair[0][1]+1, pin_second=sep_pair[0][0]+1)
        for nj, j in enumerate(sep_pair):

#            print("sep_pair= ", j, depth)
            if depth==0:
                copy_smile_backbone[sep_pair[nj][1]]=1
                copy_smile_backbone[sep_pair[nj][0]]=1
#                print("in backbone 1= ", sep_pair[nj][1], sep_pair[nj][0])
#                print("backbone inside= ", copy_smile_backbone)
            if nj!=0:
                get_smile_rings(atom_to_use,'',depth,n_rings, keep_in_mol, pin_end=sep_pair[nj][1]+1, pin_second=sep_pair[nj][0]+1)
                get_smile_rings(atom_to_use,'',depth,n_rings, keep_in_mol, pin_end=sep_pair[nj][0]+1, pin_second=sep_pair[nj][1]+1)



        sep_path=[0 for i in range(largest_sep+1)]
        sep_path[0]=sep_pair[0][0]


        if second_end!=0:
            sep_path[-1]=second_end-1
        else:
            sep_path[-1]=sep_pair[0][-1]


        for i in range(largest_sep-1):
            for j in copy_bond_partners[sep_path[i]]:
                if atom_to_use.count(j)==1 and copy_H_list[j-1]==1:
                    if copy_mol_sep_map[0][j-1][sep_path[-1]]==largest_sep-1-i:
                        sep_path[i+1]=j-1

        tempsmile+=copy_smile_type[sep_pair[0][0]]
        temp_id.append(sep_pair[0][0]+1)

        for j in copy_branch_partners[sep_path[0]]:
            if sep_path.count(j-1)==0 and atom_to_use.count(j)==1 and keep_in_mol[j-1]==1 and copy_H_list[j]==1:
                send_group=[j]+[k+1 for k in range(molsize+n_rings) if copy_mol_sep_map[0][k][j-1]<copy_mol_sep_map[0][k][sep_path[i+1]] and k!=sep_path[i+1] and k+1!=j and copy_H_list[k]==1 and keep_in_mol[k]==1]
#                print("send_group1= ", send_group)
                sub_smile=get_smile_rings(send_group,'',depth+1,n_rings, keep_in_mol, pin_end=j)
                tempsmile+="(%s)"%sub_smile[0]
                temp_id.append("(")
                temp_id.extend([k for k in sub_smile[1]])
                temp_id.append(")")

        for i in range(largest_sep):
            if depth==0:
#                print("in backbone 3= ", sep_path[i+1])
                copy_smile_backbone[sep_path[i+1]]=1
            if copy_branch_point[sep_path[i+1]]==0:
                tempsmile+=copy_smile_type[sep_path[i+1]]
                temp_id.append(sep_path[i+1]+1)

            else:
                tempsmile+=copy_smile_type[sep_path[i+1]]
                temp_id.append(sep_path[i+1]+1)
                for j in copy_branch_partners[sep_path[i+1]]:
                    if sep_path.count(j-1)==0 and atom_to_use.count(j)==1 and keep_in_mol[j-1]==1 and copy_H_list[j]==1:
                        send_group=[j]+[k+1 for k in range(molsize+n_rings) if copy_mol_sep_map[0][k][j-1]<copy_mol_sep_map[0][k][sep_path[i+1]] and k!=sep_path[i+1] and k+1!=j and copy_H_list[k]==1 and keep_in_mol[k]==1]
#                        print("send_group2= ", send_group)
                        sub_smile=get_smile_rings(send_group,'',depth+1,n_rings, keep_in_mol, pin_end=j)
                        tempsmile+="(%s)"%sub_smile[0]
                        temp_id.append("(")
                        temp_id.extend([k for k in sub_smile[1]])
                        temp_id.append(")")

    else:
        for i in atom_to_use:
            tempsmile+=copy_smile_type[i-1]
            temp_id.append(i)


    if depth==0:
        smiles.append([k for k in tempsmile])
        smiles_id.append([k for k in temp_id])
        copy_overall_smile_backbone.append(copy_smile_backbone)
#        print("we got there: ", tempsmile)
#        sys.exit()


    return [tempsmile, temp_id, depth]





















#--------------------------------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------------------------------



compress_dict={}
compress_stride={}
inv_comp_dict={}


def compress_smile_simp(smile,stride):
    print("smile= ", smile)
    comp_smile=[i for i in smile]    
    new_smile=[]
    while len(new_smile)!=len(comp_smile):
        smile_len=len(comp_smile)
        print("len of smile= ", (smile_len))
        new_smile=[i for i in comp_smile]
        for ni in range(smile_len-stride*2): 
            print("i= ", ni)
            comped=0
            comp=[comp_smile[ni+j] for j in range(stride)]
#            print("comp= ", comp)
            if contiguous_str(comp):
                print("closed comp= ", comp, smile_len, int(smile_len/stride), int((ni)/stride))
                for j in range(int(int(smile_len-ni)/stride)):
                    
                    if j>0:
                        
                        print([comp_smile[ni+k] for k in range(stride)]) 
                        print([comp_smile[ni+j*stride+k] for k in range(stride)], j)
                        if [comp_smile[ni+k] for k in range(stride)]==[comp_smile[ni+j*stride+k] for k in range(stride)]: 
                            comp.extend([comp_smile[ni+j*stride+k] for k in range(stride)])
                            comped=1
                        else:
                            break

                if comped==1:
                    temp_smile=""       
                    for j in comp_smile: temp_smile+=j
                    temp_comp=""
#                    print("total comp= ", comp)
                    for j in comp:
                        for k in j: 
                            temp_comp+=k
                    split_comp=temp_smile.split(temp_comp) 
                    if str(temp_comp) in inv_comp_dict:
                        this_comp_code=inv_comp_dict[temp_comp]
                    else:            
                        this_comp_code=comp_codes[0]
                        compress_dict.update({this_comp_code:temp_comp})
                        inv_comp_dict.update({temp_comp:this_comp_code})
                        compress_stride.update({this_comp_code:stride})
                        del comp_codes[0]

                    comp_smile=[]
                    for nj, j in enumerate(split_comp): 
                        comp_smile.extend([k for k in j])
                        if nj+1!=len(split_comp): comp_smile.append(this_comp_code)
                    break        
    return comp_smile
    


#--------------------------------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------------------------------


#smile="CAB(C)AB(C)AB(C)AB(C)AB(C)C"
#smile="CAAB(B1)BBB(B1)B(B2)BBB(B2)B(B3)BBB(B3)B(B4)BBB(B4)AAC"

#smile="DAAAAAAAAB1B(CBB1)B2CB(B(B2)AAAAAAAAD)B4CB(B(B4)AAAAAAAAD)B3BB(BS3)AAAAAAAAD"
#rev_smile="DAAAAAAAAB(BC3)BB3B(B(B4)AAAAAAAAD)CB4B(B(B2)AAAAAAAAD)CB2B(CBB1)B1AAAAAAAAD"

#smile="DEAAEAAEAB1B(CBB1)B2CB(B(B2)AEAAEAAED)B4CB(B(B4)AEAAEAAED)B3BB(BS3)AEAAEAAED"
#rev_smile="DEAAEAAEAB(BC3)BB3B(B(B4)AEAAEAAED)CB4B(B(B2)AEAAEAAED)CB2B(CBB1)B1AEAAEAAED"



def find_raster(smile,rev_smile):
    match_frac=[0 for i in smile]
    match_frac[0]=len(smile)
    smile=[i for i in smile]
    rev_smile=[i for i in rev_smile]

    ep_string_id=[]
    for ni, i in enumerate(smile):
        if ["(",")"].count(i)==0:
            ep_string_id.append(ni)
        else:
            ep_string_id.append('')   

    matching_blocks=[[] for i in range(len(smile))]
    ring_match_id=[0 for i in range(len(smile))]
    ring_match_id_rev=[0 for i in range(len(rev_smile))]

    for ni,i in enumerate(smile): 
        try:
            ring_match_id[ni]=int(i)
            smile[ni]="X"
        except ValueError:
            v=2

    for ni,i in enumerate(rev_smile):
        try:
            ring_match_id_rev[ni]=int(i)
            rev_smile[ni]="X"
        except ValueError:
            v=2

    def check_parentheses(string):
        return int( (string.count(")")-string.count("("))==0 or ((string.count(")")-string.count("("))%2==1  and (string.count("(")-string.count(")")>0)))

    def cat_string(string_vec):
        emptystr=""
        for i in string_vec: 
            emptystr+=i
        return emptystr    

    for i in range(len(smile)):
        if i!=0:
            block_match=0
            match_str=[]
            j_matches=[]
            for j in range(len(smile)):
                if j>=i:
                    if smile[j]==smile[j-i]:
                        match_str.append(smile[j])
                        j_matches.append(j)
                        match_frac[i]+=1
                    if smile[j]!=smile[j-i] or j+1==len(smile):

                        if len(match_str)-match_str.count("(")-match_str.count(")")>=1:
                            keep_string=1
                            if match_str[0]==")":        
                                del match_str[0]
                                del j_matches[0] 
                            if match_str.count("X")>=0:
                                string1_rings=[ring_match_id[k] for k in j_matches] 
                                string2_rings=[ring_match_id[k-i] for k in j_matches]
                                if len([k for k in string1_rings if k>0])>1:
                                    remove_zero_string_1=[k for k in string1_rings if k>0]
                                    remove_zero_string_2=[k for k in string2_rings if k>0]
                                    ring_dict={}
                                    for nk,k in enumerate(remove_zero_string_1):
                                        try:
                                            if ring_dict[k]!=remove_zero_string_2[nk]:
                                                keep_string=0
                                                break
                                        except KeyError:
                                            ring_dict.update({k:remove_zero_string_2[nk]})           
                            if keep_string==1: 
                                if check_parentheses(match_str)==1:
                                    if len(match_str)-match_str.count("X")-match_str.count("(")-match_str.count(")")>1:
                                        matching_blocks[i].append(match_str)
                        match_str=[]      
                        j_matches=[] 

        
        block_match=0
        match_str=[]
        j_matches=[]    
        for j in range(len(smile)):
            if 1==0:
#            if j>=i:               
                if rev_smile[j]==smile[j-i]:
                    match_str.append(rev_smile[j])
                    j_matches.append(j)

                if rev_smile[j]!=smile[j-i] or j+1==len(smile):
                    if 1==1:
                        keep_string=1
                        if match_str.count("X")>=0:
                            string1_rings=[ring_match_id_rev[k] for k in j_matches]
                            string2_rings=[ring_match_id[k-i] for k in j_matches]
                            if len([k for k in string1_rings if k>0])>1:
                                remove_zero_string_1=[k for k in string1_rings if k>0]
                                remove_zero_string_2=[k for k in string2_rings if k>0]
                                ring_dict={}
                                for nk,k in enumerate(remove_zero_string_1):
                                    try:
                                        if ring_dict[k]!=remove_zero_string_2[nk]:
                                            keep_string=0
                                            break
                                    except KeyError:
                                        ring_dict.update({k:remove_zero_string_2[nk]})
                        if keep_string==1:
                            if check_parentheses(match_str)==1:
                                if len(match_str)-match_str.count("X")-match_str.count("(")-match_str.count(")")>1:
                                    matching_blocks[i].append(match_str)
                    match_str=[]
                    j_matches=[]


        block_match=0
        match_str=[]
        j_matches=[]
        for j in range(len(smile)):
            if 1==0:
#            if j>=i:
                if rev_smile[j-i]==smile[j]:    
                    match_str.append(rev_smile[j-i])
                    j_matches.append(j)
                if rev_smile[j-i]!=smile[j] or j+1==len(smile):
                    if 1==1:
                        keep_string=1
                        if match_str.count("X")>=0:
                            string1_rings=[ring_match_id_rev[k] for k in j_matches]
                            string2_rings=[ring_match_id[k-i] for k in j_matches]
                            if len([k for k in string1_rings if k>0])>1:
                                remove_zero_string_1=[k for k in string1_rings if k>0]
                                remove_zero_string_2=[k for k in string2_rings if k>0]
                                ring_dict={}
                                for nk,k in enumerate(remove_zero_string_1):
                                    try:
                                        if ring_dict[k]!=remove_zero_string_2[nk]:
                                            keep_string=0
                                            break
                                    except KeyError:
                                        ring_dict.update({k:remove_zero_string_2[nk]})
                        if keep_string==1:
                            if check_parentheses(match_str)==1:
                                if len(match_str)-match_str.count("X")-match_str.count("(")-match_str.count(")")>1:
                                    matching_blocks[i].append(match_str)
                    match_str=[]
                    j_matches=[]

    return match_frac






def get_raster_block(vec):
    block_len=0
    if len(vec)>1:
        del vec[0]
        for ni in range(len(vec)-1):
            if ni>0:
                if vec[ni]>vec[ni-1] and vec[ni]>vec[ni+1]:
                    block_len=ni
                    break

        block_len=vec.index(max(vec))+1
        if block_len==len(vec): block_len=0
    return block_len

#--------------------------------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------------------------------


def invert_compression(smile):
    inverted=[]
    for i in smile:
        if i not in compress_dict: inverted.append([i])
        else: 
            stride=compress_stride[i]
            if stride==1:   inverted.append([compress_dict[i]])
            else:
                fully_inverted=0
                repeat=[compress_dict[i][j] for j in range(stride)]
                while fully_inverted==0:
                    fully_inverted=1
                    for nj, j in enumerate(repeat):
                        temp_j=[]
                        for nk, k in enumerate(j):
                            if k in compress_dict:
                                temp_j.extend([l for l in compress_dict[k]])
                                fully_inverted=0
                            else: temp_j.extend([k])   
                        repeat[nj]=[k for k in temp_j] 
                redone_repeat=[]
                for j in repeat: redone_repeat.extend(j) 
                inverted.append([redone_repeat for j in range(int(len(compress_dict[i])/stride))]) 
    return inverted


#--------------------------------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------------------------------
initial_schemes=[]
final_schemes=[]

def     find_scheme(n_aa, cg_group,  groups, H_list, dynamic_avg, scheme_id, n_max,*positional_parameters, **keyword_parameters):



    for i in range(10): print()

    print("cg_groups: ", cg_group)
    print("groups= ", groups)
    #Reinitialize, find which groups need help
    finished_groups=[[k for k in i] for i in groups]
#    print("n_aa, cg_group: ", n_aa, cg_group)
    mapped, terminal, groups_done = start_find(n_aa, cg_group, aa_bond_partners, finished_groups, H_list,dynamic_avg)
    print("orig stuff= ",mapped, terminal, groups_done)

    loopc=0
    while sum(mapped)!=len(cg_group):
        loopc+=1
        if loopc==20:
            print("loop break")
            sys.exit()
        groups_done=get_group_done(n_aa,cg_group,atomistic_bond_partners,finished_groups,H_list,dynamic_avg)
        print("groups_done= ", groups_done)
        #if all groups are done, add the next layer of terms 
        if sum(groups_done)==len(groups_done):
            #recalculate dynamic_avg
            for nj, j in enumerate(cg_group):
#                print("mapped= ", mapped[nj])
#                print("nj,j,terminal= ", nj, j , terminal[nj])
                if (H_list[j-1]==1 and mapped[nj]==0 and terminal[nj]==0):
                    print("adding: ", j)
                    finished_groups.append([j])
                    mapped[nj]=1



#        for nj, j in enumerate(cg_group):
#            if mapped[nj]==1:
#                for m in atomistic_bond_partners[j-1]:
#                        if cg_group.count(m)!=0:
#                            terminal[cg_group.index(m)]=0



        print("finished_groups= ", finished_groups)
        print("terminals: ", terminal)
        print("dynamic_avg= ", dynamic_avg)
#        sys.exit()

        groups_done=get_group_done(n_aa,cg_group,atomistic_bond_partners,finished_groups,H_list,dynamic_avg)
        # Find which atoms we will propose should be added
        print("groups_done= ", groups_done)
        proposed_adds=[[] for i in finished_groups]
        for ni, i in enumerate(finished_groups):
            if groups_done[ni]==0:
                for nj, j in enumerate(i):
                    for nk, k in enumerate(atomistic_bond_partners[j-1]):
                        if cg_group.count(k)!=0:
                            if H_list[k-1]==1 and mapped[cg_group.index(k)]==0 and proposed_adds[ni].count(k)==0:
                                proposed_adds[ni].append(k)
        
        if sum([len(i) for i in proposed_adds])==0:
            for nj,j in enumerate(cg_group):
                if mapped[nj]==0 and H_list[j-1]==1:
                    finished_groups.append([j])
                    mapped[nj]=1

        else:
            conflicts=[]
            conflict_beads=[]
            for nj, j in enumerate(proposed_adds):
                if len(j)>0:
                    for nk, k in enumerate(j):
                        if sum([l.count(k) for l in proposed_adds])>1 and conflict_beads.count(k)==0:
                            conflict_beads.append(k)
                            conflicts.append([l for l in range(len(finished_groups)) if proposed_adds[l].count(k)==1])

            n_config=1
            for j in conflicts: n_config=n_config*len(j)
            prefac=[prod([len(conflicts[1+j+k]) for j in range(len(conflicts)-1-k)]) for k in range(len(conflicts))]
            prefac=[prod([len(conflicts[j]) for j in range(len(conflicts))])]+prefac
            copy_of_finished=[[j for j in k] for k in finished_groups]

            n_cont=n_config
            if n_config<n_max:
                n_max=n_max/n_config
            if n_config>n_max:
                n_cont=n_max
                n_max=1

            choice_keep=random.sample([j for j in range(n_config)],n_cont)

            for nj, j in enumerate(choice_keep):
                choice_proposed=[[k for k in l] for l in proposed_adds]
                choices=[(j%prefac[k])/prefac[k+1] for k in range(len(conflicts))]
                for nc, c in enumerate(choices):
                    for nk, k in enumerate(conflicts[nc]):
                        if nk!=c:
                            choice_proposed[k].remove(conflict_beads[nc])

                n_overs=[len(k)+len(finished_groups[nk])-dynamic_avg for nk, k in enumerate(choice_proposed)]
                if sum([int(k>0) for k in n_overs])==0:
                    if nj==0:
                        for nk, k in enumerate(choice_proposed):
                            for l in k:
                                finished_groups[nk].append(l)
                                mapped[cg_group.index(l)]=1

                                for m in atomistic_bond_partners[l-1]:
                                    if cg_group.count(m)!=0:
                                        terminal[cg_group.index(m)]=0

                    else:
                        copy_of_groups=[[k for k in l] for l in copy_of_finished]
                        for nk, k in enumerate(choice_proposed):
                            for l in k:
                                copy_of_groups[nk].append(l)
                        find_scheme(n_aa, cg_group, copy_of_groups, H_list, dynamic_avg,scheme_id, n_max)

                else:
                    if nj==0:
                        # add the guys which don't go over
                        for no, o in enumerate(n_overs):
                            if groups_done[no]==0:
                                if o<=0:
                                    for l in choice_proposed[no]:
                                        finished_groups[no].append(l)
                                        mapped[cg_group.index(l)]=1
                                        for m in atomistic_bond_partners[l-1]:
                                            if cg_group.count(m)!=0:
                                                terminal[cg_group.index(m)]=0

                        # now generate the possible pairs for each over 
                        over_groups=[]
                        n_group_config=[]
                        pairs=[[] for k in range(len(finished_groups))]
                        n_add_config=1
                        for no, o in enumerate(n_overs):
                            if groups_done[no]==0:
                                if o>0:
                                    over_groups.append(no)
                                    n_group_config.append(int(math.factorial(len(choice_proposed[no]))/(math.factorial(o)*math.factorial(len(choice_proposed[no])-o))))
                                    n_add_config=n_add_config*n_group_config[-1]
                                    pairs[no]=list(itertools.permutations(choice_proposed[no],(len(choice_proposed[no])-o)))

                        over_prefac=[prod([n_group_config[1+k+l] for k in range(len(over_groups)-1-l)]) for l in range(len(over_groups))]
                        over_prefac=[prod([n_group_config[k] for k in range(len(over_groups))])]+over_prefac
                        add_copy_of_groups=[[l for l in m] for m in finished_groups]

                        n_add_max=n_max
                        n_add_cont=n_add_config

                        if n_add_config<n_add_max:
                            n_add_max=n_add_config

                        if n_add_config>n_add_max:
                            n_add_cont=n_max
                            n_add_max=1

                        choice_add_keep=random.sample([k for k in range(n_add_config)],n_add_cont)
                        for nk, k in enumerate(choice_add_keep):
                            add_choice_id=[int((k%over_prefac[l])/over_prefac[l+1]) for l in range(len(over_groups))]
                            add_choices=[pairs[over_groups[l]][add_choice_id[l]] for l in range(len(add_choice_id))]
                            if nk==0:
                                for nl, l in enumerate(add_choices):
                                    for m in l:
                                        finished_groups[over_groups[nl]].append(m)
                                        mapped[cg_group.index(m)]=1
                                        for n in atomistic_bond_partners[m-1]:
                                            if cg_group.count(n)!=0:
                                                terminal[cg_group.index(n)]=0

                            else:
                                add_2_copy_of_groups=[[l for l in m] for m in add_copy_of_groups]
                                for nl, l in enumerate(add_choices):
                                    for m in l:
                                        add_2_copy_of_groups[over_groups[nl]].append(m)
                                scheme_id+=1
                                find_scheme(n_aa, cg_group, add_2_copy_of_groups, H_list, dynamic_avg, scheme_id, n_add_max)

                    else:
                        copy_of_groups=[[k for k in l] for l in copy_of_finished]
                        for no, o in enumerate(n_overs):
                            if groups_done[no]==0:
                                if o<=0:
                                    for l in choice_proposed[no]:
                                        copy_of_groups[no].append(l)

                        over_groups=[]
                        n_group_config=[]
                        pairs=[[] for k in range(len(copy_of_groups))]
                        n_add_config=1
                        for no, o in enumerate(n_overs):
                            if groups_done[no]==0:
                                if o>0:
                                    over_groups.append(no)
                                    n_group_config.append(int(math.factorial(len(choice_proposed[no]))/(math.factorial(o)*math.factorial(len(choice_proposed[no])-o))))
                                    n_add_config=n_add_config*n_group_config[-1]
                                    pairs[no]=list(itertools.permutations(choice_proposed[no],(len(choice_proposed[no])-o)))


                        over_prefac=[prod([n_group_config[1+k+l] for k in range(len(over_groups)-1-l)]) for l in range(len(over_groups))]
                        over_prefac=[prod([n_group_config[k] for k in range(len(over_groups))])]+over_prefac
                        add_copy_of_groups=[[l for l in m] for m in copy_of_groups]

                        n_add_max=n_max
                        n_add_cont=n_add_config
                        if n_add_config<n_add_max:
                            n_add_max=n_add_config
                        if n_add_config>n_add_max:
                            n_add_cont=n_max
                            n_add_max=1
                        choice_add_keep=random.sample([k for k in range(n_add_config)],n_add_cont)

                        for nk, k in enumerate(choice_add_keep):
#                        for k in range(n_add_config):
                            add_choice_id=[int((k%over_prefac[l])/over_prefac[l+1]) for l in range(len(over_groups))]
                            add_choices=[pairs[over_groups[l]][add_choice_id[l]] for l in range(len(add_choice_id))]
                            add_2_copy_of_groups=[[l for l in m] for m in add_copy_of_groups]
                            for nl, l in enumerate(add_choices):
                                for m in l:
                                    add_2_copy_of_groups[over_groups[nl]].append(m)
                            scheme_id+=1
                            find_scheme(n_aa, cg_group, add_2_copy_of_groups, H_list, dynamic_avg, scheme_id, n_add_max)

    print("finished_add= ", finished_groups)
    initial_schemes.append(finished_groups)

#    sys.exit()
    
    return










def     find_scheme_rest(n_aa, aa_bond_partners, groups, H_list, dynamic_avg, scheme_id, n_max, *positional_parameters, **keyword_parameters):

    
#    print()
    #Reinitialize, find which groups need help
    finished_groups=[[k for k in i] for i in groups]
    mapped, terminal, groups_done = start_find(n_aa,[i+1 for i in range(n_aa)], aa_bond_partners, finished_groups, H_list,dynamic_avg)
    loopc=0

    no_merge=[0 for i in range(molsize)]
    no_merge_group=[0 for i in groups]
    

    if ('protect' in keyword_parameters): 
        no_merge=[i for i in keyword_parameters['protect']]
        no_merge_group=[int(sum([no_merge[j-1] for j in i])==len(i)) for i in groups]

#        print("groups= ", groups)
#        print("no_merge= ", no_merge)
#        print("no_mege_goup= ", no_merge_group) 
   
    while sum(mapped)!=sum(H_list):
        loopc+=1
        if loopc==20:
            print("loop break")
            sys.exit()

        groups_done=get_group_done(n_aa,[i+1 for i in range(n_aa)],atomistic_bond_partners,finished_groups,H_list,dynamic_avg)
#        print("groups= ", groups)
#        print("groups_done= ", groups_done)
        for ni,i in enumerate(no_merge_group):
            if i==1: groups_done[ni]=1
#        print("groups_done= ", groups_done)
        #if all groups are done, add the next layer of terms 
        if sum(groups_done)==len(groups_done):
            #recalculate dynamic_avg
            if len(finished_groups)<n_cg:
                v=2
            for j in range(n_aa):
                if (H_list[j]==1 and mapped[j]==0 and terminal[j]==0):
                    finished_groups.append([j+1])
                    mapped[j]=1

#        print("finished_groups= ", finished_groups)
        groups_done=get_group_done(n_aa,[i+1 for i in range(n_aa)],atomistic_bond_partners,finished_groups,H_list,dynamic_avg)
#        print("new groups_done: ", groups_done)


        for ni,i in enumerate(no_merge_group):
            if i==1: groups_done[ni]=1
#        print("groups_done= ", groups_done)

        # Find which atoms we will propose should be added
        proposed_adds=[[] for i in finished_groups]
        for ni, i in enumerate(finished_groups):
            if groups_done[ni]==0:
                for nj, j in enumerate(i):
                    for nk, k in enumerate(atomistic_bond_partners[j-1]):
                        if H_list[k-1]==1 and mapped[k-1]==0 and proposed_adds[ni].count(k)==0 and proposed_adds[ni].count(k)==0:
                            proposed_adds[ni].append(k)

#        print("proposed_adds= ", proposed_adds)

#        sys.exit()

        if sum([len(i) for i in proposed_adds])==0:
            for j in range(n_aa):
                if mapped[j]==0 and H_list[j]==1:
                    finished_groups.append([j+1])
                    mapped[j]=1
        else:
            conflicts=[]
            conflict_beads=[]
            for nj, j in enumerate(proposed_adds):
                if len(j)>0:
                    for nk, k in enumerate(j):
                        if sum([l.count(k) for l in proposed_adds])>1 and conflict_beads.count(k)==0:
                            conflict_beads.append(k)
                            conflicts.append([l for l in range(len(finished_groups)) if proposed_adds[l].count(k)==1])

#            print("conflicts", conflicts)
            n_config=1
            for j in conflicts: n_config=n_config*len(j)
            prefac=[prod([len(conflicts[1+j+k]) for j in range(len(conflicts)-1-k)]) for k in range(len(conflicts))]
            prefac=[prod([len(conflicts[j]) for j in range(len(conflicts))])]+prefac
            copy_of_finished=[[j for j in k] for k in finished_groups]

            n_cont=n_config
#            print("n_config= ", n_config)
            if n_config<n_max:
                n_max=int(n_max/n_config)
            if n_config>n_max:
                n_cont=n_max
                n_max=1

            choice_keep=random.sample([j for j in range(n_config)],n_cont)

            for nj, j in enumerate(choice_keep):
                choice_proposed=[[k for k in l] for l in proposed_adds]
                choices=[int((j%prefac[k])/prefac[k+1]) for k in range(len(conflicts))]
       
                for nc, c in enumerate(choices):
                    for nk, k in enumerate(conflicts[nc]):
                        if nk!=c:
                            choice_proposed[k].remove(conflict_beads[nc])
                n_overs=[len(k)+len(finished_groups[nk])-dynamic_avg for nk, k in enumerate(choice_proposed)]
                if sum([int(k>0) for k in n_overs])==0:
                    if nj==0:
                        for nk, k in enumerate(choice_proposed):
                            for l in k:
                                finished_groups[nk].append(l)
                                mapped[l-1]=1

                                for m in atomistic_bond_partners[l-1]:
                                    terminal[m-1]=0

                    else:
                        # generate new groups to sprout
                        copy_of_groups=[[k for k in l] for l in copy_of_finished]
                        for nk, k in enumerate(choice_proposed):
                            for l in k:
                                copy_of_groups[nk].append(l)
                        find_scheme_rest(n_aa, aa_bond_partners, copy_of_groups, H_list, dynamic_avg,scheme_id, n_max, protect=no_merge)

                else:
                    if nj==0:
                        # add the guys which don't go over
                        for no, o in enumerate(n_overs):
                            if groups_done[no]==0:
                                if o<=0:
                                    for l in choice_proposed[no]:
                                        finished_groups[no].append(l)
                                        mapped[l-1]=1
                                        for m in atomistic_bond_partners[l-1]:
                                            terminal[m-1]=0

                        # now generate the possible pairs for each over 
                        over_groups=[]
                        n_group_config=[]
                        pairs=[[] for k in range(len(finished_groups))]
                        n_add_config=1
                        for no, o in enumerate(n_overs):
                            if groups_done[no]==0:
                                if o>0:
                                    over_groups.append(no)
                                    n_group_config.append(int(math.factorial(len(choice_proposed[no]))/(math.factorial(o)*math.factorial(len(choice_proposed[no])-o))))
                                    n_add_config=n_add_config*n_group_config[-1]
                                    pairs[no]=list(itertools.permutations(choice_proposed[no],(len(choice_proposed[no])-o)))

                        over_prefac=[prod([n_group_config[1+k+l] for k in range(len(over_groups)-1-l)]) for l in range(len(over_groups))]
                        over_prefac=[prod([n_group_config[k] for k in range(len(over_groups))])]+over_prefac
                        add_copy_of_groups=[[l for l in m] for m in finished_groups]

                        n_add_max=n_max
                        n_add_cont=n_add_config

                        if n_add_config<n_add_max:
                            n_add_max=int(n_add_max/n_config)

                        if n_add_config>n_add_max:
                            n_add_cont=n_max
                            n_add_max=1

#                        print([k for k in range(n_add_config)])
                        choice_add_keep=random.sample([k for k in range(n_add_config)],n_add_cont)

                        for nk, k in enumerate(choice_add_keep):
                            add_choice_id=[int((k%over_prefac[l])/over_prefac[l+1]) for l in range(len(over_groups))]

                            add_choices=[pairs[over_groups[l]][add_choice_id[l]] for l in range(len(add_choice_id))]
                            if nk==0:
                                for nl, l in enumerate(add_choices):
                                    for m in l:
                                        finished_groups[over_groups[nl]].append(m)
                                        mapped[m-1]=1
                                        for n in atomistic_bond_partners[m-1]:
                                            terminal[n-1]=0
                            else:

                                add_2_copy_of_groups=[[l for l in m] for m in add_copy_of_groups]
                                for nl, l in enumerate(add_choices):
                                    for m in l:
                                        add_2_copy_of_groups[over_groups[nl]].append(m)
                                scheme_id+=1
                                find_scheme_rest(n_aa, aa_bond_partners, add_2_copy_of_groups, H_list, dynamic_avg, scheme_id, n_add_max, protect=no_merge)


                    else:
                        copy_of_groups=[[k for k in l] for l in copy_of_finished]
                        for no, o in enumerate(n_overs):
                            if groups_done[no]==0:
                                if o<=0:
                                    for l in choice_proposed[no]:
                                        copy_of_groups[no].append(l)

                        over_groups=[]
                        n_group_config=[]
                        pairs=[[] for k in range(len(copy_of_groups))]
                        n_add_config=1
                        for no, o in enumerate(n_overs):
                            if groups_done[no]==0:
                                if o>0:
                                    over_groups.append(no)
                                    n_group_config.append(int(math.factorial(len(choice_proposed[no]))/(math.factorial(o)*math.factorial(len(choice_proposed[no])-o))))
                                    n_add_config=n_add_config*n_group_config[-1]
                                    pairs[no]=list(itertools.permutations(choice_proposed[no],(len(choice_proposed[no])-o)))


                        over_prefac=[prod([n_group_config[1+k+l] for k in range(len(over_groups)-1-l)]) for l in range(len(over_groups))]
                        over_prefac=[prod([n_group_config[k] for k in range(len(over_groups))])]+over_prefac
                        add_copy_of_groups=[[l for l in m] for m in copy_of_groups]

                        n_add_max=n_max
                        n_add_cont=n_add_config
#                        print("orig: ", n_add_cont, n_add_max)
                        if n_add_config<n_add_max:
                            n_add_max=int(n_add_config)
                        if n_add_config>n_add_max:
                            n_add_cont=n_max
                            n_add_max=1
#                        print("n_add_config: ", n_add_config, [k for k in range(n_add_config)])
#                        print("n_add_count", n_add_cont)
                        choice_add_keep=random.sample([k for k in range(n_add_config)],n_add_cont)

                        for nk, k in enumerate(choice_add_keep):
                            add_choice_id=[int((k%over_prefac[l])/over_prefac[l+1]) for l in range(len(over_groups))]
                            add_choices=[pairs[over_groups[l]][add_choice_id[l]] for l in range(len(add_choice_id))]
                            add_2_copy_of_groups=[[l for l in m] for m in add_copy_of_groups]
                            for nl, l in enumerate(add_choices):
                                for m in l:
                                    add_2_copy_of_groups[over_groups[nl]].append(m)
                            scheme_id+=1
                            find_scheme_rest(n_aa, aa_bond_partners, add_2_copy_of_groups, H_list, dynamic_avg, scheme_id, n_add_max, protect=no_merge)
    final_schemes.append(finished_groups)
    return







def get_term_scheme(atoms_belong, dynamic_avg, is_rep):
    temp_group=[]
    if len(atoms_belong)<=dynamic_avg:
        if is_rep==1:
            temp_group=[[i for i in atoms_belong]]
        else:
            temp_group=[[[i for i in atoms_belong]]]
#        temp_group=[[[i for i in atoms_belong]]]

    else:
        find_scheme(len(atoms_belong),atoms_belong,[],H_list,dynamic_avg,69, 100)
#        print("initial_schemes= ", initial_schemes)
        temp_group=[i for i in initial_schemes]
    return temp_group




def reduce_temp(schemes, temp_belong):
    reduced_schemes=[[[i for i in j] for j in k] for k in schemes]   
 
    for ns, scheme in enumerate(reduced_schemes):

#        print(scheme)
        scheme_bond_partners=get_cg_bond_partners(scheme,aa_bond_partners,molsize)
        groups_size=[len(i) for i in scheme]
        small_size=min(groups_size)
        id_smalls=[ni for ni, i in enumerate(scheme) if len(i)==small_size]
        group_id=get_group_id(n_aa,scheme)

        small_size_nonrep=min([i for ni, i in enumerate(groups_size)])
        id_smalls_nonrep=[ni for ni, i in enumerate(scheme) if len(i)==small_size_nonrep]
        while small_size_nonrep==1:
            cur_group=random.choice(id_smalls_nonrep)
            partner_sizes=[groups_size[j-1] for j in scheme_bond_partners[cur_group]]
            smallest_neighbors=[k for nk, k in enumerate(scheme_bond_partners[cur_group]) if groups_size[k-1]==min(partner_sizes)]
            add_to=random.choice(smallest_neighbors)
#            print("adding: ", scheme[cur_group])
#            print("adding to: ", scheme[add_to-1])
            scheme[add_to-1].extend(scheme[cur_group])
            del scheme[cur_group]
            scheme_bond_partners=get_cg_bond_partners(scheme,aa_bond_partners,molsize)
            groups_size=[len(i) for i in scheme]
            small_size_nonrep=min([i for ni, i in enumerate(groups_size)])
            id_smalls_nonrep=[ni for ni, i in enumerate(scheme) if len(i)==small_size_nonrep]
#        print("reduced_scheme:", scheme)
        reduced_schemes[ns]=[[j for j in i] for i in scheme]


    def see_if_already(scheme2):
        new=1
        for k in globals()['schemes_in_progress']:
            if len(k)==len(scheme2):
                this_new=0
                for j in scheme2:
                    new2=0
                    for l in k:
                        if sorted(l)==sorted(j):
                            new2=1
                            break
                    if new2==0:
                        this_new=1
                        break
                if this_new==0:
                    new=0
                    break
        return new





    if len(reduced_schemes)>1:
        keep_schemes=[]
        for ns, scheme in enumerate(reduced_schemes):
            new=1
            for nk, k in enumerate(keep_schemes):
                if len(k)==len(scheme):
                    this_new=0
                    for nj, j in enumerate(scheme):
                        new2=0
                        for l in k:
                            if sorted(l)==sorted(j):
                                new2=1
                                break
                        if new2==0:
                            this_new=1
                            break
                    if this_new==0:
                        new=0
                        break 
            
            if new==1: keep_schemes.append([[j for j in i] for i in scheme])
        reduced_schemes=[[[j for j in i] for i in k] for k in keep_schemes]


    return reduced_schemes


#--------------------------------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------------------------------

absolute_finals=[]



rings=find_rings()[0]
print("rings= ", rings)
n_rings=len(rings)



if n_rings>0:
    print("reconfigure with ring as singular group")
    copy_bond_partners=[[j for j in i] for i in atomistic_bond_partners]

    ring_pairs=[]

    copy_H_list=[i for i in H_list]+[1 for i in rings]
    print("len H list: ", len(copy_H_list))
    copy_smile_type=[i for i in smile_type]+[0 for  i in rings]
    copy_smile_backbone.extend([0 for i in range(n_rings)])
    
    copy_branch_point=[i for i in branch_point]
    copy_branch_partners=[[j for j in i] for i in branch_partners]

    for ni, i in enumerate(rings):
        match=0
        for nj,j in enumerate(ring_pairs):
            print("test rings: ", i, ring_pairs[nj][0])
            if len(compare_aa_groups(i, ring_pairs[nj][0],atomistic_bond_partners))!=0:
                match=1 
                ring_pairs[nj].append(i)
                copy_smile_type[natom+ni]=comp_codes[nj]
                break
        if match==0:
            copy_smile_type[natom+ni]=comp_codes[len(ring_pairs)]
            ring_pairs.append([i]) 

    for i in ring_pairs: del comp_codes[0]

    print("ring_pairs= ", ring_pairs)
    print("new smile types= ", copy_smile_type)
    copy_mol_sep_map=[[[molsize+1 for i in range(molsize+n_rings)] for j in range(molsize+n_rings)]]
    keep_in_mol=[1 for i in range(natom+n_rings)]




    for ni,i in enumerate(rings):
        copy_branch_point.append(0)
        copy_bond_partners.append([])
        for j in i:
            keep_in_mol[j-1]=0
            temp_remove_vec=[]
            for k in copy_bond_partners[j-1]:
                if temp_remove_vec.count(k)==0 and i.count(k)==0:
                    temp_remove_vec.append(k)
            for nk, k in enumerate(temp_remove_vec):
                copy_bond_partners[-1].append(k)
                copy_bond_partners[k-1].append(natom+ni+1)
                copy_bond_partners[k-1].remove(j)
                copy_mol_sep_map[0][k-1][natom+ni]=1
                copy_mol_sep_map[0][natom+ni][k-1]=1
                copy_branch_partners[k-1].append(natom+ni+1)


        copy_branch_partners.append([nj+1 for nj, j in enumerate(copy_mol_sep_map[0][natom+ni]) if j==1 and copy_H_list[nj]==1])

        if len(copy_branch_partners[-1])>2:
            copy_branch_point[-1]=1
            

    for ni, i in enumerate(copy_bond_partners): 
        if keep_in_mol[ni]==1:
            copy_mol_sep_map[0][ni][ni]=0
            for nj, j in enumerate(i):
                copy_mol_sep_map[0][ni][j-1]=1


    for j in range(molsize+n_rings):
        if copy_H_list[j]==0: keep_in_mol[j]=0
        level=j+1
        for k in range(molsize+n_rings):
            if keep_in_mol[k]==1:
                katm=k
                for l in range(molsize+n_rings):
                    latm=l
                    if copy_mol_sep_map[0][k][l]==level:
                        for m in range(molsize+n_rings):
                            matm=m
                            if copy_mol_sep_map[0][l][m]==1 and m>k:
                                if level+1<copy_mol_sep_map[0][k][m]:
                                    copy_mol_sep_map[0][k][m]=level+1
                                    copy_mol_sep_map[0][m][k]=level+1


    for j in range(molsize+n_rings):
        if keep_in_mol[j]==1:
            for k in range(molsize+n_rings):
                if keep_in_mol[k]==1:
                    if j>k:
                        copy_mol_sep_map[0][j][k]=copy_mol_sep_map[0][k][j]
    

    print("going to get smiles")
    get_smile_rings([ni+1 for ni, i in enumerate( keep_in_mol) if i==1],'',0, n_rings, keep_in_mol)


    print()
    print()
    print()


    final_backbone=[0 for i in range(molsize)]

    for i in range(len(copy_overall_smile_backbone)):
#        print("overall_smile_backbone")
#        print(copy_overall_smile_backbone[i])
        for j in range(molsize): 
            if copy_overall_smile_backbone[i][j]==1:
                final_backbone[j]=1    

        for g in range(n_rings):
            if copy_overall_smile_backbone[i][natom+g]==1:
#                print("soup")
                for l in rings[g]: final_backbone[l-1]=1 
        for g in range(n_rings): del copy_overall_smile_backbone[i][-1] 
        




    final_aa_smiles=[]
    for i in smiles:
        print()
        print()
        print("aa smile= ", i)
        block_size=1
        smile_comp=compress_smile_simp(i,block_size)
        old_smile=[]
        while old_smile!=smile_comp:    
            old_smile=[j for j in smile_comp]
            rastering=find_raster(smile_comp,smile_comp)
            block_len=get_raster_block(rastering)
            if block_len>0:
                smile_comp=compress_smile_simp(smile_comp,block_len)
        final_aa_smiles.append(smile_comp)


#    print()
    print("final_aa_smiles: ")
    for i in final_aa_smiles: print(i)




    print("comp_dict= ")
    for i in compress_dict: print(i, compress_dict[i])
    print("smiles_id= ", smiles_id)







    print("final_aa_smiles= ")
    for ni, i in enumerate(final_aa_smiles):
        for j in range(10): print()
        print()
        print("INVERTING: ",i)
        inv_smile=invert_compression(i)

        in_repeat=[0 for j in range(molsize+n_rings)]
        print("inv_smile= ", inv_smile)
        # go and group the atoms
        #stores what are the "grouped" stuff


        initial_groups=[[]]
        repeating_groups=[]
        a_count=0
        for nj, j in enumerate(inv_smile):
            print()            
            list_yes=int(type(j[0]) is list)
            print("j= ", j, list_yes)
            if len(j[0])>1:
                print("this can be grouped", list_yes)
            #find what atoms belong to this
            atoms_belong=[]
            correspondence_dict=[]
            new_dynamic=dynamic_avg
            if len(initial_groups[0])<len(cg_coors[0]):
                new_dynamic=int((len(atomistic_coors[0])-a_count)/(len(cg_coors[0])-len(initial_groups[0])))
                print("new_dynamic in loop= ", new_dynamic)
            if list_yes==0:
                for k in j[0]:
                    if k!="(" and k!=")":
                        atoms_belong.append(smiles_id[ni][a_count])
                    a_count+=1

            else:
                for m in range(len(j)):
                    atoms_belong.append([])
                    for k in j[m]:
                        if k!="(" and k!=")":
                            atoms_belong[-1].append(smiles_id[ni][a_count])
                            in_repeat[smiles_id[ni][a_count]-1]=1
                        a_count+=1
                correspondence_dict=[{} for k in range(len(j))]
                for nk, k in enumerate(atoms_belong[0]):
                    if k>molsize:
                        # this is gonna be a ring correspondence
                        for l in range(len(j)):
                            this_corres=compare_aa_groups(rings[k-natom-1],rings[atoms_belong[l][nk]-natom-1])[0]                           
                            for m in this_corres:
                                correspondence_dict[l].update({m: this_corres[m]})
    
                    else:
                        for l in range(len(j)):
                            correspondence_dict[l].update({k:atoms_belong[l][nk]})

                repeating_groups.append(correspondence_dict)



            for nk, k in enumerate(rings):
                if in_repeat[molsize+nk]==1:
                    for nl, l in enumerate(k):
                        in_repeat[l-1]=1 




            is_rep=0
            if len(correspondence_dict)>0:
                temp_belong=[k for k in atoms_belong[0]]
                is_rep=1
            else:
                temp_belong=[k for k in atoms_belong]

            temp_ring_replace=[k for k in temp_belong if k>molsize]
            for k in temp_ring_replace:
                temp_belong.remove(k)
                temp_belong.extend(rings[k-molsize-1])

            # END ATOM BELONG BLOCK



            initial_schemes=[]


            temp_scheme=get_term_scheme(temp_belong, dynamic_avg, is_rep)
            copy_initial_groups=[[[k for k in l] for l in m] for m in initial_groups]



            print("j= ", j)
            if len(j[0])>1:
                print("old_temp: ", temp_scheme, is_rep)
                
#                temp_scheme=reduce_temp(temp_scheme,temp_belong)
                print("new_temp: ", temp_scheme)
                new_initial=[]
                for nk, k in enumerate(temp_scheme):
                    print("temp= ", k)
                    copy2_initial_groups=[[[n for n in l] for l in m] for m in copy_initial_groups]
                    temp_list_yes=int(type(k[0]) is list)
                    print("temp_list_yes: ", temp_list_yes)
                    if list_yes==0:
                        for l in copy2_initial_groups:
                            if temp_list_yes==1:
                                for n in k:
                                    l.append([m for m in n])
                            else:
                                l.append([m for m in k])
                    else:


                        if temp_list_yes==0:
#                        print("corresponding add: ", k)
                            for nn, n in enumerate(copy2_initial_groups):
#                            for nx, x in enumerate(k):
#                                print("x= ", x)
                                for nl, l in enumerate(correspondence_dict):
                                    n.append([correspondence_dict[nl][m] for m in k])
                        else:
                            print("going in")
                            for nn, n in enumerate(copy2_initial_groups):
                                print("n= ", n)
                                for nl, l in enumerate(correspondence_dict):
                                    print("l= ", l)
                                    for nm, m in enumerate(k):
                                        print("m= ", m, [x for x in m])
                                        n.append([l[x] for x in m])
#                                    print("CORRES: ", [correspondence_dict[nl][m] for m in x])
                    for l in copy2_initial_groups:
                        new_initial.append([[m for m in n] for n in l])
                initial_groups=[[[m for m in l] for l in k] for k in new_initial]
             
#            print("initial_gorups= ")
#            for k in initial_groups: print(k) 

        final_schemes=[]
#        print("initial_gorups= ")
#        for j in initial_groups: 
#            print()
#            print(j)

        for j in initial_groups:
#            print("j= ",j)
#            print("in_repeat= ", in_repeat)
            find_scheme_rest(n_aa, aa_bond_partners, j, H_list, dynamic_avg, 69, 100, protect=in_repeat)
#            print(final_schemes[0])
    
        print()
        print()        
#        print("overall smiles= ")
#        print(copy_overall_smile_backbone)
        for nj, j in enumerate(final_schemes):
            print()
            print("nj= ", nj)
            print(j)
            for l in repeating_groups: print(l)


            suffix="final_scheme_%s_%s"%(ni,nj)
            make_colored_trj(n_aa, atomistic_coors, H_list, j,lens,suffix)
            fileo.fileo_trans([[nl+1]+[k for k in l] for nl,l in enumerate(j)],["" for l in range(len(j)+1)],"initial_scheme_%s_%s.txt"%(ni+1,nj+1))
            fileo.fileo_trans([[nl+1]+[k for k in l] for nl,l in enumerate(repeating_groups)],["repeating_groups"],"initial_scheme_%s_%s.txt"%(ni+1,nj+1), append=True)
            fileo.fileo([[k+1 for k in range(len(final_backbone))],final_backbone],["id","in_backbone"],"initial_scheme_%s_%s.txt"%(ni+1,nj+1), append=True)
































































else:
    get_smile([i+1 for i in range(molsize)],'',0)


    print("smiles= ")
    for i in smiles: print(i)


#    sys.exit()

     
    final_aa_smiles=[]
    final_cg_smiles=[]

    for i in smiles:
#        for j in range(10): print()
#        print("og smile= ", i)
         
        block_size=1
        print("block_size= ", block_size)
        smile_comp=compress_smile_simp(i,block_size)
        print("smile_comp= ", smile_comp)
        old_smile=[]
        while old_smile!=smile_comp:
            old_smile=[j for j in smile_comp]
            rastering=find_raster(smile_comp,smile_comp)
            block_len=get_raster_block(rastering)   
            print("block_len= ", block_len)

            if block_len>0:
                smile_comp=compress_smile_simp(smile_comp,block_len)       
            print("comped smile= ", smile_comp)

        final_aa_smiles.append(smile_comp)



    for i in cg_smiles:
        block_len=1
        print("block_len= ", block_len)
        smile_comp=compress_smile_simp(i,block_len)
        old_smile=[]
        while old_smile!=smile_comp:
            old_smile=[j for j in smile_comp]
            rastering=find_raster(smile_comp,smile_comp)
            block_len=get_raster_block(rastering)
            print("block_len= ", block_len)
            if block_len>0:
                smile_comp=compress_smile_simp(smile_comp,block_len)
        final_cg_smiles.append(smile_comp)


    final_backbone=[0 for i in range(molsize)]

    for i in range(len(overall_smile_backbone)):
        for j in range(molsize):
            if overall_smile_backbone[i][j]==1:
                final_backbone[j]=1



#    print("smile_id", smiles_id)
#    for i in smiles_id: print(i)
#    print("final_aa_smiles= ")
#    for ni, i in enumerate(final_aa_smiles):
#        print("i= ", i)


#    print("comp_dict= ", compress_dict)

    for ni, i in enumerate(final_aa_smiles): 
#        for j in range(10):print()
        print("i= ", i)

        inv_smile=invert_compression(i)
#        print("inv_smile= ", inv_smile)

        in_repeat=[0 for j in range(molsize)]
        
        
        # go and group the atoms
        #stores what are the "grouped" stuff
        initial_groups=[[]]
        repeating_groups=[]
        a_count=0

        smiles_id_copy=[j for j in smiles_id[ni]]
        print("smiles_id_copy= ", smiles_id_copy)
        repeat_units=[j for j in inv_smile if int(type(j[0]) is list)==1]
        non_repeat_units=[j for j in inv_smile if int(type(j[0]) is list)!=1]
        units=[repeat_units,non_repeat_units]
        print("repeat_units= ", repeat_units)
        print("non_repeats= ", non_repeat_units)
#        sys.exit()
    
        n_mapped=0

#        sys.exit()      
 
#        for nu, u in enumerate(units):             
#            for nj, j in enumerate(u):
        if 1==1:
            for nj, j in enumerate(inv_smile):
                list_yes=int(type(j[0]) is list)
                print("j= ", j)
                if len(j[0])>1:
                    print("this can be grouped", list_yes)
                #find what atoms belong to this
                atoms_belong=[]
                correspondence_dict=[]
                new_dynamic=dynamic_avg
                if len(initial_groups[0])<len(cg_coors[0]): 
                    new_dynamic=int((len(atomistic_coors[0])-a_count)/(len(cg_coors[0])-len(initial_groups[0])))
                print("new_dynamic sup= ", new_dynamic)
                print("smiles_id= ", smiles_id)
#                dynamic_avg=new_dynamic
                if list_yes==0:
                    for k in j[0]:
                        if k!="(" and k!=")":
                            atoms_belong.append(smiles_id[ni][a_count])
                        a_count+=1 
                    
                else:
                    for m in range(len(j)):
                        atoms_belong.append([])
                        for k in j[m]:
                            if k!="(" and k!=")":
                                atoms_belong[-1].append(smiles_id[ni][a_count])
                                in_repeat[smiles_id[ni][a_count]-1]=1
                            a_count+=1
                    correspondence_dict=[{} for k in range(len(j))]
                    for nk, k in enumerate(atoms_belong[0]):   
                        for l in range(len(j)):
                            correspondence_dict[l].update({k:atoms_belong[l][nk]}) 
    #                print("correspondence_dict= ", correspondence_dict)
                    repeating_groups.append(correspondence_dict)


            
                is_rep=0
                if len(correspondence_dict)>0:
                    temp_belong=[k for k in atoms_belong[0]]     
                    is_rep=1
                else:
                    temp_belong=[k for k in atoms_belong]

                initial_schemes=[]        
       
                 
                # calculated dynamic_avg
                print("temp_belong= ", temp_belong)
                print("dynamic_avg= ", dynamic_avg) 
                temp_scheme=get_term_scheme(temp_belong, dynamic_avg, is_rep)
                print("temp_scheme= ", temp_scheme)

#                sys.exit()

                copy_initial_groups=[[[k for k in l] for l in m] for m in initial_groups]
                

                if len(j[0])>1:            
                    new_initial=[]    
                    print("temp_scheme= ", temp_scheme)
                    for nk, k in enumerate(temp_scheme):
                        print()
                        copy2_initial_groups=[[[n for n in l] for l in m] for m in copy_initial_groups]
                        print("copy2_initial: ", copy2_initial_groups)
                        print("k= ", k)
                        temp_list_yes=int(type(k[0]) is list)
    #                    print("temp list yes: ", temp_list_yes, list_yes)
                        if list_yes==0:
                            for l in copy2_initial_groups:
                                if temp_list_yes==1:
                                    for n in k:
                                        l.append([m for m in n])      
                                else: 
                                    l.append([m for m in k])
                        else:
                            if temp_list_yes==0:
                                for n in copy2_initial_groups:
                                    for nl, l in enumerate(correspondence_dict):
                                        n.append([correspondence_dict[nl][m] for m in k]) 
                            else:
                                for n in copy2_initial_groups:
                                    for nl, l in enumerate(correspondence_dict):
                                        for nm, m in enumerate(k):
                                            n.append([l[x] for x in m])


                    for k in copy2_initial_groups: 
                        new_initial.append(k)
                    initial_groups=[[[m for m in l] for l in k] for k in new_initial]
                

                
                print("initial_groups= ", initial_groups)

                
                
        final_schemes=[]
        for j in initial_groups: 
            find_scheme_rest(n_aa, aa_bond_partners, j, H_list, dynamic_avg, 69, 100, protect=in_repeat)


        for nj, j in enumerate(final_schemes):
            print()
            print("nj= ", nj)
            print(j)
            for l in repeating_groups: print(l)

            suffix="final_scheme_%s_%s"%(ni,nj)        
            make_colored_trj(n_aa, atomistic_coors, H_list, j,lens,suffix)
            fileo.fileo_trans([[nl+1]+[k for k in l] for nl,l in enumerate(j)],["" for l in range(len(j)+1)],"initial_scheme_%s_%s.txt"%(ni+1,nj+1))
            fileo.fileo_trans([[nl+1]+[k for k in l] for nl,l in enumerate(repeating_groups)],["repeating_groups"],"initial_scheme_%s_%s.txt"%(ni+1,nj+1), append=True)        

            print([k+1 for k in range(len(final_backbone))])
            print(final_backbone)

            fileo.fileo([[k+1 for k in range(len(final_backbone))],final_backbone],["id","in_backbone"],"initial_scheme_%s_%s.txt"%(ni+1,nj+1), append=True)

#        sys.exit()





