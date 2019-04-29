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

# Constants
gram2Dalton = 6.02213665e23
m2Ang = 1.0e10
cm2Ang = 1.0e8

# State Data
state={}
state['N']=0.0
state['M']=0.0
state['rhoma']=0.0
state['rhogcm']=0.0
state['a']=1.0
state['b']=1.0
state['c']=1.0

def print_head():
	print("""______                    _  _            _____         _               _         _                
|  _  \                  (_)| |          /  __ \       | |             | |       | |               
| | | |  ___  _ __   ___  _ | |_  _   _  | /  \/  __ _ | |  ___  _   _ | |  __ _ | |_   ___   _ __ 
| | | | / _ \| '_ \ / __|| || __|| | | | | |     / _` || | / __|| | | || | / _` || __| / _ \ | '__|
| |/ / |  __/| | | |\__ \| || |_ | |_| | | \__/\| (_| || || (__ | |_| || || (_| || |_ | (_) || |   
|___/   \___||_| |_||___/|_| \__| \__, |  \____/ \__,_||_| \___| \__,_||_| \__,_| \__| \___/ |_|   
                                   __/ |                                                           
                                  |___/ """)

def print_help():
	print_head()
	print("""
A small utility for doing density-related calculations for MM simulations.

Command line arguments:
	(none)                  - By default, DC will enter interactive mode.
	--lammpsdata file.data  - Reads the box dimensions and atoms in the LAMMPS data file
	                          and returns the density in different units.
	""")

def densfromLAMMPS(fn):
	f=open(fn,'r')
	natoms=-1
	ntypes=-1
	atomsDone=False
	massesDone=False
	xread=False
	yread=False
	zread=False
	masses={}
	atoms=[]
	nmols=0
	#read the header
	lines=0
	while(True):
		l=f.readline().lower()
		lines+=1
		if('atoms' in l):
			natoms=int(l.split()[0])
		elif('atom types') in l:
			ntypes=int(l.split()[0])
		elif('xlo') in l:
			ls=l.split()
			x=np.abs(float(ls[1])-float(ls[0]))
			xread=True
		elif('ylo') in l:
			ls=l.split()
			y=np.abs(float(ls[1])-float(ls[0]))
			yread=True
		elif('zlo') in l:
			ls=l.split()
			z=np.abs(float(ls[1])-float(ls[0]))
			zread=True
		if((natoms>0)and(ntypes>0)and(xread)and(yread)and(zread)):
			break
	print("We have %d atoms from %d types."%(natoms,ntypes))
	print("The volume of the simulation box is %0.3f angs**3"%(x*y*z))
	print("We read %d lines thus far."%(lines))
	#read the body
	while(True):
		l=f.readline()
		if(l.strip().lower()=='masses'):
			l=f.readline() #empty line
			for n in range(ntypes):
				l=f.readline().split()
				masses[l[0]]=float(l[1])
			massesDone=True
		elif((len(l.strip())>=5)and(l.strip().lower()[:5]=='atoms')):
			l=f.readline() #empty line
			for n in range(natoms):
				l=f.readline().split()
				atoms.append(l[2])
				if(int(l[1])>nmols):
					nmols=int(l[1])
			atomsDone=True
		if((massesDone)and(atomsDone)):
			break #go on to processing the data
	f.close()
	print("We read %d atoms from %d types, forming %d molecules."%(len(atoms),len(masses.keys()),nmols))
	mass=0.0
	for a in atoms:
		mass += masses[a]
	print("Total mass of the system is %0.3f dalton (%0.3e g)"%(mass,mass/gram2Dalton))
	print("The density of the box is:")
	print(" %15.5f molecules.ang**-3"%(nmols/(x*y*z)))
	print(" %15.5f dalton.ang**-3"%(mass/(x*y*z)))
	print(" %15.5f g.cm**-3"%((mass/gram2Dalton)/((x/cm2Ang)*(y/cm2Ang)*(z/cm2Ang))))
	print(" %15.5f kg.m**-3"%((mass/(1000*gram2Dalton))/((x/m2Ang)*(y/m2Ang)*(z/m2Ang))))

def mainMenu():
	while(True):
		print("0 - Quit")
		print("1 - Calculate density given M, N and box size.")
		opt=input("> ")
		try:
			o=int(opt)
		except:
			print("Invalid Option: %s"%(opt))
			continue
		if(o in list(range(2))):
			return(o)
		else:
			print("Invalid Option: %s"%(opt))

def rhoFromMNSize():
	while(True):
		print("Current State:")
		print("Molecular Mass: %12.4f Da"%(state['M']))
		print("Number of Molecules: %d"%(state['N']))
		print("Box Dimensions: %6.2f x %6.2f x %6.2f angs."%(state['a'],state['b'],state['c']))
		print("---")
		nmols=state['N']
		mass=state['M']*nmols
		x=state['a']
		y=state['b']
		z=state['c']
		print("Total mass of the system is %0.3f dalton (%0.3e g)"%(mass,mass/gram2Dalton))
		print("The density of the box is:")
		print(" %15.5f molecules.ang**-3"%(nmols/(x*y*z)))
		print(" %15.5f dalton.ang**-3"%(mass/(x*y*z)))
		print(" %15.5f g.cm**-3"%((mass/gram2Dalton)/((x/cm2Ang)*(y/cm2Ang)*(z/cm2Ang))))
		print(" %15.5f kg.m**-3"%((mass/(1000*gram2Dalton))/((x/m2Ang)*(y/m2Ang)*(z/m2Ang))))
		print("---")
		print("0 - Quit to main menu.")
		print("1 - Change Molecular Mass")
		print("2 - Change Number of Molecules")
		print("3 - Change Box Dimensions")
		opt=input(">> ")
		try:
			o=int(opt)
		except:
			print("Invalid Option: %s"%(opt))
			continue
		if(o==0):
			return 0
		elif(o==1):
			state['M']=float(input("Insert new molecular mass (in Da) >> "))
			continue
		elif(o==2):
			state['N']=int(input("Insert number of molecules >> "))
			continue
		elif(o==3):
			opt=input("Box dimensions (3 numbers, or one number for cubic box) >> ")
			opt=opt.split()
			if(len(opt)==3):
				state['a']=float(opt[0])
				state['b']=float(opt[1])
				state['c']=float(opt[2])
			else:
				state['a']=float(opt[0])
				state['b']=float(opt[0])
				state['c']=float(opt[0])
			continue


def mainLoop():
	print_head()
	while(True):
		opt=mainMenu()
		if(opt==0):
			break
		elif(opt==1):
			rhoFromMNSize()

if(__name__=='__main__'):
	if(len(sys.argv)<2):
		#enter interactive mode
		mainLoop()
		sys.exit(0)
	if(('-h' in sys.argv)or('--help' in sys.argv)):
		print_help()
		sys.exit(0)
	#Read arguments for non-interactive operation
	n=1
	while(n<len(sys.argv)):
		arg=sys.argv[n]
		if(arg=='--lammpsdata'):
			n +=1
			fn=sys.argv[n]
			print("Calculating density for simulation box in %s..."%(fn))
			densfromLAMMPS(fn)
			sys.exit(0)
		else:
			print("Unknown argument: %s"%(arg))
			sys.exit(1)


