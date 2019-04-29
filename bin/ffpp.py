#! /usr/bin/env python3
# -*- coding: utf8 -*-

##########################################################################
#                                                                        #
# Programa:                                 Data: __/__/____             #
#                                                                        #
# Uso:                                                                   #
#                                                                        #
# _DESCRICAO_                                                            #
#                                                                        #
#                                                                        #
#                                                                        #
#                                                                        #
# (c) Filipe Teixeira                                                    #
#                                                                        #
##########################################################################

import sys
import numpy as np

Bohr2Ang=0.529

Symbols=['h','he','li','be','b','c','n','o','f','ne',
'na','mg','al','si','p','s','cl','ar', 'k','ca','sc',
'ti','v','cr','mn','fe','co','ni','cu','zn','ga','ge',
'as','se','br','kr', 'rb', 'sr', 'y', 'zr', 'nb', 'mo', 
'tc','ru', 'rh', 'pd', 'ag', 'cd', 'in', 'sn', 'sb', 'te',
'i', 'xe', 'cs', 'ba', 'la', 'ce', 'pr', 'nd', 'pm', 'sm', 
'eu', 'gd', 'tb', 'dy', 'ho', 'er', 'tm', 'yb', 'lu',
'hf', 'ta', 'w', 're', 'os', 'ir', 'pt', 'au', 'hg', 'tl', 'pb', 'bi']

Masses=[1.007825, 4.002602, 6.94, 9.0121831, 10.81, 12.0000, 14.007, 15.9949159, 18.998403163, 20.1797, 22.98976928, 24.305, 26.9815385, 28.085, 30.973761998, 32.06, 35.45, 39.948, 39.0983, 40.078, 44.955908, 47.867, 50.9415, 51.9961, 54.938044, 55.845, 58.933194, 58.6934, 63.546, 65.38, 69.723, 72.63, 74.921595, 78.971, 79.904, 83.798, 85.4678, 87.62, 88.90584, 91.224, 92.90637, 95.95, 97, 101.07, 102.9055, 106.42, 107.8682, 112.414, 114.818, 118.71, 121.76, 127.6, 126.90447, 131.293, 132.90545196, 137.327, 138.90547, 140.116, 140.90766, 144.242, 145, 150.36, 151.964, 157.25, 158.92535, 162.5, 164.93033, 167.259, 168.93422, 173.054, 174.9668, 178.49, 180.94788, 183.84, 186.207, 190.23, 192.217, 195.084, 196.966569, 200.592, 204.38, 207.2, 208.9804]

Raddi=[ 0.5, 0.3, 1.7, 1.1, 0.9, 0.7, 0.6, 0.6, 0.5, 0.4, 1.9, 1.5, 1.2, 1.1, 1.0, 0.9, 0.8, 0.7, 2.4, 1.9, 1.8, 1.8, 1.7, 1.7, 1.6, 1.6, 1.5, 1.5, 1.5, 1.4, 1.4, 1.3, 1.1, 1.0, 0.9, 0.9, 2.7, 2.2, 2.1, 2.1, 2.0, 1.9, 1.8, 1.8, 1.7, 1.7, 1.7, 1.6, 1.6, 1.5, 1.3, 1.2, 1.2, 1.1, 3.0, 2.5, 2.0, 2.0, 2.5, 2.1, 2.1, 2.4, 2.3, 2.3, 2.3, 2.3, 2.3, 2.3, 2.2, 2.2, 2.2, 2.1, 2.0, 1.9, 1.9, 1.9, 1.8, 1.8, 1.7, 1.7, 1.6, 1.5, 1.4]


## Input/Output ##
def readXYZplain(ifn): 
	""" Reads and returns the data in a XYZ file
	Modified to also read labels, charges and vdw parameters from a modified
	xyz file, in which each line has the following format:
	Symbol XX YY ZZ [Label] [Charge][vdw Epsilon] [vdw Sigma] 
	Where the last two parameters are given in Angstrons and kcal/mol, respectively"""
	f=open(ifn,'r')
	data=f.readlines()
	f.close()
	natoms=int(data[0])
	s=[]
	g=np.zeros((natoms,3))
	labels=[]
	charges=np.zeros(natoms)
	sigma=np.zeros(natoms)
	epsilon=np.zeros(natoms)
	i=0
	for n in range(2,2+natoms):
		l=data[n].split()
		if(len(l)>=5):
			labels.append(l[4])
		if(len(l)>=6):
			charges[i]=float(l[5])
		if(len(l)==8):
			sigma[i]=float(l[7])
			epsilon[i]=float(l[6])
		s.append(l[0])
		g[i,:]=np.array(list(map(float,l[1:4])))
		i += 1
	return(s,g,labels,charges,sigma,epsilon)

def readXYZ(ifn): 
	""" Reads and returns the data in a XYZ file
	Modified to also read labels, charges and vdw parameters from a modified
	xyz file, in which each line has the following format:
	Symbol XX YY ZZ [Label] [Charge] [vdw Epsilon] [vdw Sigma] 
	Where the last two parameters are given in kcal/mol and Angstrons, respectively.
	Outputs a dic"""
	o={}
	f=open(ifn,'r')
	data=f.readlines()
	f.close()
	natoms=int(data[0])
	s=[]
	g=np.zeros((natoms,3))
	labels=[]
	hasCharges=False
	hasVDW=False
	charges=np.zeros(natoms)
	sigma=np.zeros(natoms)
	epsilon=np.zeros(natoms)
	i=0
	for n in range(2,2+natoms):
		l=data[n].split()
		if(len(l)>=5):
			labels.append(l[4])
		if(len(l)>=6):
			charges[i]=float(l[5])
			hasCharges=True
		if(len(l)==8):
			epsilon[i]=float(l[6])
			sigma[i]=float(l[7])
			hasVDW=True
		s.append(l[0])
		g[i,:]=np.array(list(map(float,l[1:4])))
		i += 1
	o['symbols']=s
	o['geometry']=g
	if(labels!=[]):
		o['labels']=labels
	if(hasCharges):
		o['charges']=charges
	if(hasVDW):
		o['sigma']=sigma
		o['epsilon']=epsilon
	return(o)

def readLabels(ifn):
	f=open(ifn,'r')
	data=f.readlines()
	f.close()
	natoms=int(data[0])
	labs=[]
	for n in range(2,2+natoms):
		labs.append(data[n].split()[-1].strip()) # We'll allow possible numbering of the lines :)
	return(labs)

def readFloats(ifn):
	f=open(ifn,'r')
	data=f.readlines()
	f.close()
	natoms=int(data[0])
	o=np.zeros(natoms)
	i=0
	for n in range(2,2+natoms):
		# We'll allow possible numbering of the lines :)
		o[i]=float(data[n].split()[-1].strip())
		i+=1
	return(o)


## Support functions ##
def makeLabels(ll):
	o=[]
	fmt="%s%02d"
	if(len(ll)>100):
		fmt="%s%03d"
	for n in range(len(ll)):
		o.append(fmt%(ll[n].capitalize(),n+1))
	return(o)

def bondLength(geo,bl):
	"""Returns distance between atoms bl[1] and bl[0]
	(numbering starts at 0) in geometry geo"""
	return np.linalg.norm(geo[bl[1]]-geo[bl[0]])

def angleAmp(geo,al,deg=False):
	"""Returns amplitude (optnialy in degs) for the valence angle formed by
	the atoms in al, with al[1] being the apex (numbering starts at 0)."""
	r01=geo[al[0]]-geo[al[1]]
	r21=geo[al[2]]-geo[al[1]]
	r01=r01/np.linalg.norm(r01)
	r21=r21/np.linalg.norm(r21)
	phi=np.arccos(np.dot(r01,r21))
	if(deg):
		phi=np.rad2deg(phi)
	return phi

def oopAmp(geo,al,deg=False):
	"""Returns amplitude (optnialy in degs) for the out-of-plane angle formed by
	the atoms in al, with al[1] being the central atom (numbering starts at 0)."""
	r1=geo[al[0]]-geo[al[1]]
	r2=geo[al[2]]-geo[al[1]]
	r3=geo[al[3]]-geo[al[1]]
	r1 /= np.linalg.norm(r1)
	r2 /= np.linalg.norm(r2)
	r3 /= np.linalg.norm(r3)
	n1=np.cross(r1,r2)
	n2=np.cross(r3,r2)
	atmp=np.dot(n1,n2)/(np.linalg.norm(n1)*np.linalg.norm(n2))
	if(atmp>1.0):
		phi=np.arccos(1.0)
	elif(atmp<-1.0):
		phi=np.arccos(-1.0)
	else:
		phi=np.arccos(atmp)
	# r1 is normal to the reference plane, so
	# project n1 and n2 onto the reference plane
	n1p = n1 - (np.dot(n1,r1)*r1)
	n2p = n2 - (np.dot(n2,r1)*r1)
	nc=np.cross(n1p,n2p)
	if(np.dot(nc,r1)>0.0):
		phi = (2.0*np.pi)-phi
	if(deg):
		phi=np.rad2deg(phi)
	return phi

def torsionAmp(geo,al,deg=False):
	"""Returns amplitude (optnialy in degs) for the 0-1-2-3 torsion angle formed by
	the atoms in al, with al[1] and al[2] defining the central bond
	(numbering starts at 0)."""
	r1=geo[al[0]]-geo[al[1]]
	r2=geo[al[2]]-geo[al[1]]
	r3=geo[al[3]]-geo[al[2]]
	r1 /= np.linalg.norm(r1)
	r2 /= np.linalg.norm(r2)
	r3 /= np.linalg.norm(r3)
	n1=np.cross(r1,r2)
	n2=np.cross(r3,r2)
	if(np.linalg.norm(n1)<1.0e-3):
		# possible hazard for linear bonds!!!
		return 0.0
	elif(np.linalg.norm(n2)<1.0e-3):
		return 0.0
	n1 /= np.linalg.norm(n1)
	n2 /= np.linalg.norm(n2)
	atmp=np.dot(n1,n2)/(np.linalg.norm(n1)*np.linalg.norm(n2))
	if(atmp>1.0):
		phi=np.arccos(1.0)
	elif(atmp<-1.0):
		phi=np.arccos(-1.0)
	else:
		phi=np.arccos(atmp)
	# r2 is normal to the reference plane, so
	# project on n1 and n2 onto the reference plane
	n1p = n1 - (np.dot(n1,r2)*r2)
	n2p = n2 - (np.dot(n2,r2)*r2)
	nc=np.cross(n1p,n2p)
	if(np.dot(nc,r1)>0.0):
		phi = (2.0*np.pi)-phi
	if(deg):
		phi=np.rad2deg(phi)
	return phi

def makeBondList(molecule,tol=0.3):
	"""Autoatic determination of the connectivity between atoms"""
	o=[] # Each bond is a tupple with the two atomic indexes
	symbols=molecule['symbols']
	geo=molecule['geometry']
	natoms=len(symbols)
	bondList=[]
	for i in range(natoms):
		r1=Raddi[Symbols.index(symbols[i].lower())]
		for j in range(i+1,natoms):
			rl=(1.0+tol)*(r1+Raddi[Symbols.index(symbols[j].lower())])
			if((bondLength(geo,[i,j])<=rl)and((i,j) not in bondList)):
				bondList.append((i,j))
	return(bondList)

def makeConnectivity(system, bondList):
	"""Returns a list of the connections for each atom, ordered by relevance"""
	o=[]
	s=system['symbols']
	natoms=len(s)
	#first make a list of lists with the connectivity for each atom
	for i in range(natoms):
		o.append([])
	for i in range(natoms):
		for j in range(i+1,natoms):
			if((i,j) in bondList):
				o[i].append(j)
				o[j].append(i)
	# now, for each atom, order the list
	ol=[]
	for i in range(natoms):
		ul=[(Masses[Symbols.index(s[j].lower())],j) for j in o[i]]
		ul.sort()
		ul.reverse()
		ol.append([i[1] for i in ul])
	return(ol)

def appendIC(icl,c):
	#print(icl,c)
	if(c[0]==1): # bond
		c1=[c[0],c[1],c[2],c[3],c[4],c[5]]
		c2=[c[0],c[1],c[3],c[2],c[4],c[5]]
		if((c1 not in icl)and(c2 not in icl)):
			icl.append(c)
		return(icl)
	elif(c[0]==2): # angle
		c1=[c[0],c[1],c[2],c[3],c[4],c[5]]
		c2=[c[0],c[1],c[4],c[3],c[2],c[5]]
		if((c1 not in icl)and(c2 not in icl)):
			icl.append(c)
		return(icl)
	elif(c[0]==3): # improper
		c1=[c[0],c[1],c[2],c[3],c[4],c[5]]
		c2=[c[0],c[1],c[2],c[3],c[5],c[4]]
		c3=[c[0],c[1],c[5],c[3],c[4],c[2]]
		c4=[c[0],c[1],c[4],c[3],c[2],c[5]]
		c5=[c[0],c[1],c[5],c[3],c[2],c[4]]
		c6=[c[0],c[1],c[4],c[3],c[5],c[2]]
		if((c1 not in icl)and(c2 not in icl)and(c3 not in icl)and
		(c4 not in icl)and(c5 not in icl)and(c6 not in icl)):
			icl.append(c)
		return(icl)
	elif(c[0]==4): # dihedral
		c1=[c[0],c[1],c[2],c[3],c[4],c[5]]
		c2=[c[0],c[1],c[5],c[4],c[3],c[2]]
		if((c1 not in icl)and(c2 not in icl)):
			icl.append(c)
		return(icl)

def makeIClist(system,connectivity):
	"""Creates a list of Internal coordinates based on the
	connectivity of the atoms, aiming at a 3N-6 objective"""
	o=[] # each IC is a tupple of 6 ints: a class, a type, follwed by the atomic indexes
	# Classes are:
	# 1 - Bonds
	# 2 - Angles
	# 3 - Impropers
	# 4 - Dihedrals
	# types depend on the atoms involved and are attributed elsewhere
	natoms=len(system['symbols'])
	for i in range(natoms):
		if(len(connectivity[i])==0):
			print("ERROR: Unconnected atom %d!"%(i+1))
			sys.exit(1)
		elif(len(connectivity[i])==1):
			#Ignore terminal atoms
			continue
		elif(len(connectivity[i])==2):
			#centre of a 3-atom cluster: 2 bonds + 1 angle
			o=appendIC(o,[1,0,i,connectivity[i][0],0,0])
			o=appendIC(o,[1,0,i,connectivity[i][1],0,0])
			o=appendIC(o,[2,0,connectivity[i][0],i,connectivity[i][1],0])
		elif(len(connectivity[i])==3):
			#centre of a 4-atom cluster: 3 bonds + 2 angles + 1 improper
			o=appendIC(o,[1,0,i,connectivity[i][0],0,0])
			o=appendIC(o,[1,0,i,connectivity[i][1],0,0])
			o=appendIC(o,[1,0,i,connectivity[i][2],0,0])
			o=appendIC(o,[2,0,connectivity[i][0],i,connectivity[i][1],0])
			o=appendIC(o,[2,0,connectivity[i][0],i,connectivity[i][2],0])
			o=appendIC(o,[3,0,connectivity[i][0],i,connectivity[i][1],connectivity[i][2]])
		elif(len(connectivity[i])==4):
			#centre of a 5-atom cluster: check if this is an isolated molecule.
			o=appendIC(o,[1,0,i,connectivity[i][0],0,0])
			o=appendIC(o,[1,0,i,connectivity[i][1],0,0])
			o=appendIC(o,[1,0,i,connectivity[i][2],0,0])
			o=appendIC(o,[1,0,i,connectivity[i][3],0,0])
			o=appendIC(o,[2,0,connectivity[i][0],i,connectivity[i][1],0])
			o=appendIC(o,[2,0,connectivity[i][0],i,connectivity[i][2],0])
			o=appendIC(o,[2,0,connectivity[i][0],i,connectivity[i][3],0])
			o=appendIC(o,[2,0,connectivity[i][1],i,connectivity[i][2],0])
			allTerm=True
			for j in connectivity[i]:
				if(len(connectivity[j])==1):
					allTerm=False
			if(allTerm):
				#single 5-atom molecule: 4 bonds + 4 angles + 1 improper
				o=appendIC(o,[3,0,connectivity[i][0],i,connectivity[i][1],connectivity[i][2]])
			else:
				#5-atom custer centred on i: 4 bonds + 4 angles + 1 dihedral
				#find a reference atom from the connectivity of the dihedral
				ref=-1
				for j in connectivity[i]:
					for k in connectivity[j]:
						if(k!=i):
							ref=k
							for l in connectivity[i]:
								if(l not in [i,j,k]):
									o=appendIC(o,[4,0,k,j,i,l])
									break
						if(ref!=-1):
							break
					if(ref!=-1):
						break
		else:
			print("Warnning: connectivity for n=%d not yet implemented!"%(len(connectivity[i])))
			o.append([])
	return(o)

def labelsFromConnects(s,c):
	"""Generates labels from connectivity"""
	natoms=len(s['symbols'])
	o=[]
	for i in range(natoms):
		con=len(c[i])
		symb=s['symbols'][i].capitalize()
		sl=[symb]
		if(symb=='H'): # H's are labeled after the attoms they connect to
			o.append("H")
		else:
			o.append("%s_%dval"%(symb.capitalize(),con))
	# now do the H's
	for i in range(natoms):
		symb=s['symbols'][i].capitalize()
		if(symb=='H'): # H's are labeled after the attoms they connect to
			o[i]="H%s"%(o[c[i][0]])
	return(o)

def classifyIC(s,icl):
	"""Determines the type of each internal coordinate based on the
	labels for each atom. """
	icnames=[]
	types=[[],[],[],[]] 
	l=s['labels']
	for i in range(len(icl)):
		icType="UNK"
		icClass=icl[i][0]
		if(icClass==1):
			icType="B_%s-%s"%(l[icl[i][2]],l[icl[i][3]])
			if(icType not in types[0]):
				types[0].append(icType)
			icl[i][1]=types[0].index(icType)+1
		elif(icClass==2):
			icType="A_%s-%s-%s"%(l[icl[i][2]],l[icl[i][3]],l[icl[i][4]])
			if(icType not in types[1]):
				types[1].append(icType)
			icl[i][1]=types[1].index(icType)+1
		elif(icClass==3):
			icType="I_%s-%s-%s-%s"%(l[icl[i][2]],l[icl[i][3]],l[icl[i][4]],l[icl[i][5]])
			if(icType not in types[2]):
				types[2].append(icType)
			icl[i][1]=types[2].index(icType)+1
		elif(icClass==4):
			icType="D_%s-%s-%s-%s"%(l[icl[i][2]],l[icl[i][3]],l[icl[i][4]],l[icl[i][5]])
			if(icType not in types[3]):
				types[3].append(icType)
			icl[i][1]=types[3].index(icType)+1
		icnames.append(icType) #name the coordinate
	return(icl,icnames)

def output(fn,s,icl,icns):
	"""Outputs systems characteristics to an output pre-processsing file"""
	f=open(fn,'w')
	f.write("System Built by ffpp 1.0, the ForceField Pre-Processor\n")
	f.write("Written by Filipe Teixeira\n")
	f.write("\nThis software is provided as-is whitout any warranty, under the GPL2 license.\n\n")
	# Atomic data
	natoms=len(s['symbols'])
	f.write("\nKEY:Atoms:%d\n"%(natoms))
	f.write("# ID Type Symbol X Y Z Mass Charge Epsilon Sigma Label \n")
	atypes=[]
	for i in range(len(s['symbols'])):
		atomType=-1
		if(s['labels'][i] not in atypes):
			atypes.append(s['labels'][i])
		atomType = atypes.index(s['labels'][i])+1
		line="%4s %4d %10.3f %10.3f %10.3f %8.2f "%(s['symbols'][i], atomType, s['geometry'][i,0],s['geometry'][i,1],s['geometry'][i,2],Masses[Symbols.index(s['symbols'][i].lower())])
		if('charges' in s.keys()):
			line += "%+8.3f "%(s['charges'][i])
		else:
			line += "%8s "%('UNK')
		if('epsilon' in s.keys()):
			line += "%+8.3f "%(s['epsilon'][i])
		else:
			line += "%8s "%('UNK')
		if('sigma' in s.keys()):
			line += "%+8.3f "%(s['sigma'][i])
		else:
			line += "%8s "%('UNK')
		line += "%10s\n"%(s['labels'][i])
		f.write(line)
	# Internal coordinates
	n=np.zeros(4)
	for i in range(len(icl)):
		n[icl[i][0]-1]+=1
	# Bonds
	f.write("\nKEY:Bonds:%d\n"%(n[0]))
	f.write("# ID Type I1 I2 Kr r0 Label\n")
	count=0
	for i in range(len(icl)):
		if(icl[i][0]==1):
			count += 1
			f.write("%4d %4d %4d %4d %10s %8.3f %20s\n"%(count,icl[i][1],icl[i][2],icl[i][3],
			'UNK',bondLength(s['geometry'],icl[i][2:4]),icns[i][2:]))
	# Angles
	f.write("\nKEY:Angles:%d\n"%(n[1]))
	f.write("# ID Type I1 I2 I3 K theta0 Label\n")
	count=0
	for i in range(len(icl)):
		if(icl[i][0]==2):
			count += 1
			f.write("%4d %4d %4d %4d %4d %10s %8.1f %20s\n"%(count,icl[i][1],icl[i][2],
			icl[i][3], icl[i][4], 'UNK',angleAmp(s['geometry'],icl[i][2:5],True),icns[i][2:]))
	# Impropers
	f.write("\nKEY:Impropers:%d\n"%(n[2]))
	f.write("# ID Type I1 I2 I3 I4 K Chi0 Label\n")
	count=0
	for i in range(len(icl)):
		if(icl[i][0]==3):
			count += 1
			f.write("%4d %4d %4d %4d %4d %4d %10s %8.1f %20s\n"%(count,icl[i][1],icl[i][2],
			icl[i][3], icl[i][4], icl[i][5], 'UNK',oopAmp(s['geometry'],icl[i][2:6],True),
			icns[i][2:]))
	# Dihedrals
	f.write("\nKEY:Dihedrals:%d\n"%(n[3]))
	f.write("# ID Type I1 I2 I3 I4 K1 K2 K3 K4 Phi0 Label\n")
	count=0
	for i in range(len(icl)):
		if(icl[i][0]==4):
			count += 1
			f.write("%4d %4d %4d %4d %4d %4d %5s %5s %5s %5s %7.1f %20s\n"%(count,icl[i][1],
			icl[i][2], icl[i][3], icl[i][4], icl[i][5], 'UNK','UNK','UNK','UNK',
			torsionAmp(s['geometry'],icl[i][2:6],True), icns[i][2:]))
	f.close()

if(__name__=='__main__'):
	#handle options TODO
	if(len(sys.argv)<2):
		print("""   ___    ___    ___    ___ 
  / __\  / __\  / _ \  / _ \\
 / _\   / _\   / /_)/ / /_)/
/ /    / /    / ___/ / ___/ 
\/     \/     \/     \/     
                            
Usage: %s file(s).xyz

Where each input file generates a different forcefield pre-processed file, ffpp.

FFPP can accept ordinary XYZ files, or a slightly modified version with a header 
on the first line containning the number of atoms, an ignored comment line and 
then one line per atom with the following fields:

SYMBOL X Y Z LABEL CHARGE VDW_EPSILON VDW_SIGMA

where all distances (X, Y, Z and VDW_SIGMA) are given in angstrongs, charges 
are given in atomic units and energies in kcal/mol.
 """%(sys.argv[0]))
	for fn in sys.argv[1:]:
		#read labels, charges, etc
		molecule=readXYZ(fn)
		#make bond list
		bondList=makeBondList(molecule)
		#determine connectivity for each atom
		#Order the neighbourhood of each atom by priority
		connects=makeConnectivity(molecule,bondList)
		#Determine the ICs defined around each non-terminal atom
		icl=makeIClist(molecule,connects)
		icl.sort()
		if('labels' not in molecule.keys()):
			molecule['labels']=labelsFromConnects(molecule,connects)
		icl,icnames=classifyIC(molecule,icl)
		#Generate files
		output(fn[:-4]+'.ffpp',molecule,icl,icnames)

