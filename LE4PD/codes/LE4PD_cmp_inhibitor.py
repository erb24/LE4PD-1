# coding: utf-8

def gen_protinfo(PROTNAME,G96,TOP):
	import subprocess
	import numpy as np

	protname = str(PROTNAME)
	N = int(subprocess.check_output('grep -c "CA" '+str(TOP),shell=True))
	NFRS = int(subprocess.check_output('grep -c "TIMESTEP" '+str(G96),shell=True))
	NATOMS = int(subprocess.check_output('grep -c "ATOM" '+str(TOP),shell=True))

	array = np.array([protname,N,NFRS,NATOMS],dtype=str)
	np.savetxt("protname.txt",array.T,fmt="%s")

def get_chain_info(PROTNAME,G96, TOP):
	import subprocess
	import numpy as np
	from string import ascii_uppercase

	subprocess.call("grep '^ATOM' " + TOP +" > tmp",shell=True)

	#Going to use a kludge to perform the chain assignment, since there
	#are really nine versus six chains in the system.

	#Start with chain "A"
	counter = 0
	chain_idx = ascii_uppercase[counter]
	data = []
	num = 0
	with open("tmp","r+") as f:
		for line in f:
			if num == 0:
				dummy = int(line[22:26])
			#Skip the CA atoms in the ligand
			if line.split()[3] == 'lig_L':
				pass
			#Skip the CA atoms in the methylated C-terminal capping group
			elif line.split()[3] == 'NMA':
				pass
			else:
				#Check for 'jumps' in the sequence, indicating a new chain
				if (abs(int(line[22:26]) - dummy) >= 2) and num != 0:
					print(int(line[22:26]))
					counter += 1
					chain_idx = ascii_uppercase[counter]
				data.append(line.replace(line[20:23], line[20] + chain_idx + line[22]))
			dummy = int(line[22:26])
			num += 1


	chain_list = []
	atype_list = []
	for i in range(len(data)):
		chain_list.append(data[i][21])
		atype_list.append(data[i].split()[2])
	chain_list = np.array(chain_list)

	sum_list = []
	for c in ascii_uppercase:
		sum_list.append((chain_list == str(c)).sum())
		
	sum_list = np.array(sum_list)

	NATOMS = 0
	for i in range(len(sum_list)):
		if sum_list[i] != 0:
			#np.savetxt('natoms'+str(i+1)+'.dat',np.array([sum_list[i]]))
			with open('natoms'+str(i+1)+'.dat','w+') as f:
				f.write(str(int(sum_list[i]))+"\n")
				NATOMS += int(sum_list[i])
			f.close()
			
	nres_list = np.zeros(len(ascii_uppercase))
	for i,chain in enumerate(chain_list):
		if atype_list[i] == "CA":
			for num,c in enumerate(ascii_uppercase):
				if chain == str(c):
					nres_list[num] += 1
					continue
					
	NRES = 0
	for i in range(len(nres_list)):
		if sum_list[i] != 0:
			#np.savetxt('nres'+str(i+1)+'.dat',np.array([nres_list[i]]))
			with open('nres'+str(i+1)+'.dat','w+') as f:
				f.write(str(int(nres_list[i]))+"\n")
				NRES += int(nres_list[i])
			f.close()

	with open('nmol.dat','w+') as f:
		f.write(str((sum_list != 0).sum())+"\n")
	f.close()
	#np.savetxt('nmol.dat',np.array([(sum_list != 0).sum()],dtype=int))
	subprocess.call("rm -rfv tmp",shell=True)

	#Fix the order of the topology and trajectory
	resnum_list = []
	for num,i in enumerate(data):
		if i.split()[2] == "CA": resnum_list.append(int(i[22:26]))

	resnum_list = np.array(resnum_list)


	#TODO: Need to change the hard-coded 1147 here since the inhibitor-bound simulations 
	#use only the receptor-binding domain.

	#EDIT: For some reason, these simulations use structures with the correct ordering in the PDB file.
	#Comment out instead.

	# In[6]:

	'''res_counter = nres_list[np.where(nres_list != 0)[0]].cumsum()
	where_list = np.where(resnum_list > 1147)[0]


	# In[7]:

	shuffle_list = []
	dummy_list = list(where_list)
	head_list = resnum_list[np.where(resnum_list <= 1147)]
	counter = 0
	shuffle_list.append(counter)
	counter += 1
	hist = True
	while counter < len(resnum_list):
	#for i in range(1,len(head_list)):
		#if dummy_list == []:
		#    shuffle_list.append(counter)
		#    counter += 1
		#    hist = True
		#    continue
		if (resnum_list[counter] == resnum_list[counter - 1] + 2) and (resnum_list[counter] < 1147) and (hist == True):
			#data_array = np.insert(dummy, i, )
			shuffle_list.append(dummy_list[0])
			dummy_list.pop(0)
			hist = False
		elif (resnum_list[counter] > 1147):
			counter += 1 
			continue
		else:
			shuffle_list.append(counter)
			counter += 1
			hist = True
	'''

	# In[8]:

	CA_list = []
	for i in range(len(data)):
		if data[i].split()[2] == "CA": CA_list.append(data[i])


	# In[9]:

	with open('CA.pdb','w+') as f:
		for i in np.array(CA_list): #[shuffle_list]:
			f.write(i+'\n')
	f.close()


	# In[10]:

	#txt = np.genfromtxt('protname.txt',dtype=str)
	#protname = txt[0]
	#N = int(txt[1])
	#nfrs = int(txt[2])
	#natoms = int(txt[3])
	nmol = int(np.loadtxt('nmol.dat'))

	with open("protname.txt","w+") as f:
		f.write(PROTNAME + "\n")
		f.write(str(NRES) + "\n")
		NFRS = int(subprocess.check_output('grep -c "TIMESTEP" '+str(G96),shell=True))
		f.write(str(NFRS) + "\n")
		f.write(str(NATOMS) + "\n")
	f.close()
	print(PROTNAME, NRES, NFRS, NATOMS ,nmol)

	#Allocate arrays containing the number of alpha-carbons per molecule
	nca = np.zeros(nmol)
	for i in range(nmol):
			nca[i] = int(np.loadtxt('nres'+str(i+1)+".dat"))
	nca_cumsum = np.cumsum(nca)

	#print(protname,N,nfrs,natoms,nmol)

	N = NRES
	nfrs = NFRS

	rx = np.zeros((N,nfrs))
	ry = np.zeros((N,nfrs))
	rz = np.zeros((N,nfrs))
	lx = np.zeros((N - nmol, nfrs))
	ly = np.zeros((N - nmol, nfrs))
	lz = np.zeros((N - nmol, nfrs))
	Rinv = np.zeros((N,N))
	traj = np.load('unformatted_traj.npy')
	for numba,k in enumerate(range(0,N*nfrs,N)):
			rx[:,numba] = traj[k:k+N,0]
			ry[:,numba] = traj[k:k+N,1]
			rz[:,numba] = traj[k:k+N,2]

	#rx = rx[shuffle_list]
	#ry = ry[shuffle_list]
	#rz = rz[shuffle_list]


	# In[11]:

	#Define bond vectors
	counter = 0
	mol_counter = 0
	i = 1
	while i <= N:
		   if i == nca_cumsum[mol_counter]:
				   i += 1
				   mol_counter += 1
				   if mol_counter == nmol: break
		   else:
				   lx[counter,:] = rx[i,:] - rx[i-1,:]
				   ly[counter,:] = ry[i,:] - ry[i-1,:]
				   lz[counter,:] = rz[i,:] - rz[i-1,:]
				   counter += 1
				   i += 1


	# In[12]:

	lavm = np.zeros(N - nmol)
	lavmsq = np.zeros(N - nmol)
	avdot = np.zeros(N - nmol)
	avblsq = 0
	avbl = 0

	avgx = lx.mean(1)
	avgy = ly.mean(1)
	avgz = lz.mean(1)

	for i in range(N - nmol):
		for k in range(nfrs):
			dummy = lx[i,k]**2 + ly[i,k]**2 + lz[i,k]**2
			lavm[i] += np.sqrt(dummy)
			lavmsq[i] += dummy
			avblsq += dummy
		avbl += lavm[i]
		avdot[i] = avgx[i]**2 + avgy[i]**2 + avgz[i]**2
	lavm = lavm/nfrs
	lavmsq = lavmsq/nfrs
	avblsq = (avblsq/(N - nmol)) / nfrs
	avbl = (avbl/(N - nmol)) / nfrs

	print(avbl,avblsq)

	shuffled_traj = np.column_stack([rx.ravel(order='F'),ry.ravel(order='F'),rz.ravel(order='F')])

	subprocess.call("cp -v unformatted_traj.npy unshuffled_traj.npy", shell=True)
	np.save('unformatted_traj.npy',shuffled_traj)

'''def gen_protinfo(PROTNAME,PDB,TOP):
		import subprocess
		import numpy as np

		protname = str(PROTNAME)
		N = int(subprocess.check_output('grep -c "CA" '+str(TOP),shell=True))
		NFRS = int(subprocess.check_output('grep -c "MODEL" '+str(PDB),shell=True)-1)
		NATOMS = int(subprocess.check_output('grep -c "ATOM" '+str(TOP),shell=True))

		array = np.array([protname,N,NFRS,NATOMS],dtype=str)
		np.savetxt("protname.txt",array.T,fmt="%s")'''

'''def convert_traj(traj):
	import numpy as np
	traj = np.loadtxt(str(traj))
	#Save the trajectory in units of nm
	np.save('unformatted_traj.npy',traj/10.0)'''

def convert_traj(G96):
	import numpy as np
	import subprocess
	subprocess.call("sed '/BOX/, +1 d' "+str(G96)+" | sed '/TITLE/, +1 d' | awk 'NF==3' > tmp",shell=True)
	traj = np.loadtxt('tmp')

	np.save('unformatted_traj.npy',traj)
	subprocess.call('rm -rfv tmp',shell=True)

def make_unformatted_traj(G96):
	import numpy as np
	import sys

	txt = np.genfromtxt('protname.txt',dtype=str)
	protname = txt[0]
	N = int(txt[1])
	nfrs = int(txt[2])
	#natoms = int(txt[3])

	print(G96)
	print(protname,N,nfrs)

	dummy = True                    
	linr = []
	with open(G96) as f:
		for line in f:
			value = line.split()[0].isalpha()
			if (dummy == False) and (value == False):
				linr.append(line_prev.split())
				linr.append(line.split())
				for i in range(N-2):
					linr.append(f.readline().split())
			dummy = value
			line_prev = line

	np.save('unformatted_traj.npy',np.array(linr,dtype=float))

def Umatrix(traj):
	import numpy as np
	txt = np.genfromtxt('protname.txt',dtype=str)
	protname = txt[0]
	N = int(txt[1])
	nfrs = int(txt[2])
	natoms = int(txt[3])
	nmol = int(np.loadtxt('nmol.dat'))

	#Allocate arrays containing the number of alpha-carbons per molecule
	nca = np.zeros(nmol)
	for i in range(nmol):
		nca[i] = int(np.loadtxt('nres'+str(i+1)+".dat"))
	nca_cumsum = np.cumsum(nca)


	print(protname,N,nfrs,natoms)

	rx = np.zeros((N,nfrs))
	ry = np.zeros((N,nfrs))
	rz = np.zeros((N,nfrs))
	lx = np.zeros((N - nmol, nfrs))
	ly = np.zeros((N - nmol, nfrs))
	lz = np.zeros((N - nmol, nfrs))
	Rinv = np.zeros((N,N))
	traj = np.load(str(traj))
	for numba,k in enumerate(range(0,N*nfrs,N)):
		rx[:,numba] = traj[k:k+N,0]
		ry[:,numba] = traj[k:k+N,1]
		rz[:,numba] = traj[k:k+N,2]

	#Define bond vectors
	counter = 0
	mol_counter = 0
	i = 1
	while i <= N:
		if i == nca_cumsum[mol_counter]:
			i += 1
			mol_counter += 1
			if mol_counter == nmol: break
		else:
			lx[counter,:] = rx[i,:] - rx[i-1,:] 
			ly[counter,:] = ry[i,:] - ry[i-1,:]
			lz[counter,:] = rz[i,:] - rz[i-1,:]
			counter += 1
			i += 1

	#Calculate average inverse distances for the
	#hydrodynamic interaction matrix
	for i in range(N):
			for j in range(i,N):
				if i == j:
					Rinv[i,i] = np.nan
				else:
					Rinv[i,j] = (1/np.sqrt((rx[i,:] - rx[j,:])**2 + (ry[i,:] - ry[j,:])**2 + (rz[i,:] - rz[j,:])**2)).sum()
				Rinv[j,i] = Rinv[i,j]
	Rinv = Rinv/nfrs
	lavm = np.zeros(N - nmol)
	lavmsq = np.zeros(N - nmol)
	avdot = np.zeros(N - nmol)
	avblsq = 0
	avbl = 0

	avgx = lx.mean(1)
	avgy = ly.mean(1)
	avgz = lz.mean(1)

	for i in range(N - nmol):
		for k in range(nfrs):
			dummy = lx[i,k]**2 + ly[i,k]**2 + lz[i,k]**2
			lavm[i] += np.sqrt(dummy)
			lavmsq[i] += dummy
			avblsq += dummy 
		avbl += lavm[i]
		avdot[i] = avgx[i]**2 + avgy[i]**2 + avgz[i]**2
	lavm = lavm/nfrs
	lavmsq = lavmsq/nfrs
	avblsq = (avblsq/(N - nmol))/nfrs
	avbl = (avbl/(N - nmol))/nfrs

	print(avbl,avblsq)

	Umat = np.zeros((N - nmol,N - nmol))
	for i in range(N - nmol):
		for j in range(N - nmol):
			Umat[i,j] = (np.dot(lx[i,:],lx[j,:]) + np.dot(ly[i,:],ly[j,:]) + np.dot(lz[i,:],lz[j,:]))/(lavm[i]*lavm[j]*nfrs)

	np.save('Umatrix.npy',Umat)
	np.savetxt('Umatrix',np.insert(np.ravel(Umat),N - nmol,0))
	np.save('Rinv.npy',Rinv)
	np.savetxt('Rij',np.ravel(Rinv))
	np.savetxt('length',lavm)
	np.savetxt('lengthsq',lavmsq)
	np.savetxt('avldot.dat',avdot)
	np.savetxt('avbl',np.array([avbl]))
	np.savetxt('avblsq',np.array([avblsq]))

	#return avbl,avblsq,Umat

def fric_calc(PROTNAME,TOP):
	import numpy as np
	import sys
	import os
	import subprocess

	#TOP = str(sys.argv[1])
	pi = np.pi


	#Get basic information from the protname.txt file
	txt = np.genfromtxt('protname.txt',dtype=str)
	protname = txt[0]
	N = int(txt[1])
	nfrs = int(txt[2])
	natoms = int(txt[3])

	print(protname,N,nfrs,natoms)

	#Calculate the Miller radius per bead
	mradlist = []
	#with open(TOP) as f:
	with open ("CA.pdb") as f:
		for line in f:
			if line[0:4] != 'ATOM':
				#print(line)
				pass
			elif line[0:4] == 'ATOM' and line.split()[2] == "CA":
				dummy = line.split()

				#Really horrendous and vestigial; probably smoother
				#to make a dictionary with the residue names plus their
				#assoicated Miller radii
				if dummy[3] == "ALA": mradlist.append((113.0/(4*pi))**.5)
				elif dummy[3] == "ARG" : mradlist.append((241.0/(4*pi))**.5)
				elif dummy[3] == "ASN" : mradlist.append((158.0/(4*pi))**.5)
				elif dummy[3] == "ASP" : mradlist.append((151.0/(4*pi))**.5)
				elif dummy[3] == "CYS" : mradlist.append((140.0/(4*pi))**.5)
				elif dummy[3] == "GLN" : mradlist.append((189.0/(4*pi))**.5)
				elif dummy[3] == "GLU" : mradlist.append((113.0/(4*pi))**.5)
				elif dummy[3] == "GLY" : mradlist.append((85.0/(4*pi))**.5)
				elif dummy[3] == "HIS" : mradlist.append((194.0/(4*pi))**.5)
				elif dummy[3] == "ILE" : mradlist.append((182.0/(4*pi))**.5)
				elif dummy[3] == "LEU" : mradlist.append((180.0/(4*pi))**.5)
				elif dummy[3] == "LYS" : mradlist.append((211.0/(4*pi))**.5)
				elif dummy[3] == "MET" : mradlist.append((204.0/(4*pi))**.5)
				elif dummy[3] == "PHE" : mradlist.append((218.0/(4*pi))**.5)
				elif dummy[3] == "PRO" : mradlist.append((143.0/(4*pi))**.5)
				elif dummy[3] == "SER" : mradlist.append((122.0/(4*pi))**.5)
				elif dummy[3] == "THR" : mradlist.append((146.0/(4*pi))**.5)
				elif dummy[3] == "TRP" : mradlist.append((259.0/(4*pi))**.5)
				elif dummy[3] == "TYR" : mradlist.append((229.0/(4*pi))**.5)
				elif dummy[3] == "VAL" : mradlist.append((160.0/(4*pi))**.5)
				elif dummy[3] == "NLM" : mradlist.append((158.0/(4*pi))**.5) #ASN-NAG
				elif dummy[3] == "NMG" : mradlist.append((158.0/(4*pi))**.5) #ASN-NAG
				elif dummy[3] == "CYX" : mradlist.append((140.0/(4*pi))**.5) #Some modified CYS

	array = np.array(mradlist,dtype=str)
	np.savetxt('mrad.dat',np.array(mradlist).T,fmt='%s')


	#Calculate the average solvent-exposed surface area per bead
	if os.path.exists("resarea.xvg"):
		pass
	else:
		subprocess.call("echo '1' | gmx_mpi sasa -f "+str(PROTNAME)+".xtc -s "+str(TOP)+" -or resarea.xvg -dt 100",shell=True)
	resarea = []

	#Ignore 'residue' 26 (ACE on N-terminal resiude) and the second of 'residue' 1147 (NME on C-terminal residue)
	chk = False
	with open('resarea.xvg') as f:
		for line in f:
			if line[0] == '#' or line[0] == '@':
				pass
			else:
				#if line.split()[0] == '26':
				#	continue
				#if line.split()[0] == dummy:
				#	continue
				resarea.append(float(line.split()[1]))
			dummy = line[0]

	#Re-order the residues in the output to account for the non-canonical residues
	#Or not
	#shuffle_list = np.loadtxt('shuffle_list.dat')

	#darray = np.zeros_like(shuffle_list,dtype=int)
	#for i, numba in enumerate(shuffle_list):
	#	darray[i] = int(shuffle_list[i])
	#shuffle_list = darray
	#print(darray)
	#print(shuffle_list)

	#resarea = np.array(resarea)[shuffle_list]
	#resarea = list(resarea)

	rad = []
	for area in resarea:
		rad.append(((area/(4*np.pi))**0.5)*10)

	np.savetxt('avresrad',np.array(rad),fmt="%f")
	fratio = (np.array(rad).sum()/N)/10
	print('fratio: ',fratio)

	np.savetxt('fratio',np.array([fratio]))

	#Calculate the friction coefficients

	kB = 1.38066E-23
	try:
		T = float(np.loadtxt('temp'))
		print('Temperature (K): ',T)
	except OSError:
		print('Temperature not set. Defaulting to 300 K')
		T = 300
		np.savetxt('temp',np.array([T]))

	try:
		intv = float(np.loadtxt('internalv'))
		print('Internal viscosity factor: ',intv)
	except OSError:
		print('Internal viscosity not set. Defaulting to 2.71828')
		intv = 2.71828
		np.savetxt('intv',np.array([intv]))

	#Use NIST formula for viscosity -- good NEAR room temperature and physiological.
	#Won't work higher than, say, 360 K.

	try:
		#Load viscosity and convert to Pa s
		viscosity = float(np.loadtxt('visc.txt'))
	except OSError:
		print('No viscosity given. Using the NIST formula, which is only valid for physiological conditions,')
		print('i.e. between about 273 and 310 K.')
		viscosity = (.2131590-1.96290E-3*T+(.00246411*T)**2+(-.0018462*T)**3)
		np.savetxt('visc.txt',np.array([viscosity/1000]))

	print("Viscosity (Pa s): ",viscosity)
	

	try:
		fd20 = float(np.loadtxt('fd20'))
	except OSError:
		fd20 = 0
		np.savetxt('fd20',np.array([fd20]))

	rv = np.array(mradlist)
	rw = np.array(rad)
	rp = np.zeros(N)
	friw = np.zeros(N)
	fri = np.zeros(N)
	friwt = 0
	frit = 0
	for i in range(N):
		if rw[i] < rv[i]: 
			rp[i] = (rv[i]**2 - rw[i]**2)**0.5
		else:
			rp[i] = 0

		friw[i] = 6.0*pi*(rw[i]/10)*viscosity
		fri[i] = 6.0*pi*(rp[i]/10)*(intv*viscosity) + 6.0*pi*(rw[i]/10)*viscosity
		friwt += friw[i]
		frit += fri[i]

	avfr = frit/float(N)
	avfrw = friwt/float(N)
	np.savetxt('avfr',np.array([avfr*1.0E-9]))

	avblsq = float(np.loadtxt('avblsq'))
	sigma = (3*kB*T*1E15)/(avblsq*avfr)

	with open('sigma','w') as f:
		f.write('sigma, 1/ps\n')
		f.write(str(sigma)+'\n')

	with open('sigma.dat','w') as f:
		f.write(str(sigma)+'\n')

	fric = np.zeros((N+1,2))

	fric[0,0] = avfrw
	fric[0,1] = avfr
	for i in range(N):
		fric[i+1,:] = np.column_stack([friw[i],fri[i]])

	np.savetxt('fric',fric)

	return fratio,sigma,fric

def LUI_calc(fratio,avblsq,sigma,fric,Rinv):
	import numpy as np
	#Get basic information from the protname.txt file
	txt = np.genfromtxt('protname.txt',dtype=str)
	protname = txt[0]
	N = int(txt[1])
	nfrs = int(txt[2])
	natoms = int(txt[3])
	nmol = int(np.loadtxt('nmol.dat'))

	#Allocate arrays containing the number of alpha-carbons per molecule
	nca = np.zeros(nmol)
	for i in range(nmol):
		nca[i] = int(np.loadtxt('nres'+str(i+1)+".dat"))
	nca_cumsum = np.cumsum(nca)

	print(protname,N,nfrs,natoms)

	avfr = fric[0,1]

	M = np.zeros((N,N))

	idx = 0
	counter = 1
	for i in range(N):
		if i == idx:
			for j in range(N):
				M[i,j] = 1/N
			idx = nca_cumsum[counter - 1]
			counter += 1
		else:
			for j in range(N):
				if i == j + 1:
					M[i,j] = -1
				elif i == j:
					M[i,j] = 1

	a = np.zeros((N - nmol,N))
	counter = 0
	for i in range(N):
		if ((i == 0) or (np.where(nca_cumsum == i)[0].size == 1)):
			#print('here')
			pass
		else:
			a[counter,:] = M[i,:]
			counter += 1

	H = np.zeros((N,N))

	for i in range(N):
		for j in range(i,N):
			if i == j:
				H[i,i] = avfr/fric[i+1,1]
			else:
				H[i,j] = fratio*Rinv[i,j]
				H[j,i] = H[i,j]

	L = np.matmul(a,np.matmul(H,a.T))

	LINV = np.linalg.inv(L)
	L = np.matmul(a,np.matmul(H,a.T))
	LINV = np.linalg.inv(L)
	UINV = np.load("Umatrix.npy")
	UILI = np.matmul(UINV,LINV)

	#Eigendecomposition of UILI to find eigenvalues and eigenvectors 
	eigval,Q = np.linalg.eig(UILI)
	eigval = 1/eigval
	QINV = np.linalg.inv(Q)

	#Implement fratio loop to avoid negative eigenvalues of LU matrix
	while (np.any(1/eigval <= 1e-9) and fratio >= 0.0):
		fratio -= 0.01
		print(fratio,(1/eigval)[np.where(1/eigval <= 0)])
		H = np.zeros((N,N))
		for i in range(N):
			for j in range(i,N):
				if i == j:
					H[i,i] = avfr/fric[i+1,1]
				else:
					H[i,j] = fratio*Rinv[i,j]
					H[j,i] = H[i,j]

		L = np.matmul(a,np.matmul(H,a.T))

		LINV = np.linalg.inv(L)
		L = np.matmul(a,np.matmul(H,a.T))
		LINV = np.linalg.inv(L)
		UINV = np.load("Umatrix.npy")
		UILI = np.matmul(UINV,LINV)

		#Eigendecomposition of UILI to find eigenvalues and eigenvectors 
		eigval,Q = np.linalg.eig(UILI)
		eigval = 1/eigval
		QINV = np.linalg.inv(Q)

	perm = np.argsort(np.abs(eigval))

	Q_sorted = np.copy(Q)[:,perm]
	QINV_sorted = np.linalg.inv(Q_sorted)

	eigval_sorted = abs(eigval)[perm]

	mu = 1/(np.diag(np.matmul(QINV_sorted,np.matmul(UINV,QINV_sorted.T))))

	np.save("UILImatrix",UILI)
	np.savetxt('UILImatrix',np.ravel(UILI))
	np.save('Lmatrix',L)
	np.save("Hmatrix",H)
	np.save("Qmatrix.npy",Q_sorted)
	np.savetxt("Qmatrix",np.ravel(Q_sorted))
	np.save("QINVmatrix.npy",QINV_sorted)
	np.savetxt("QINVmatrix",np.ravel(QINV_sorted))
	np.save("lambda_eig.npy",eigval_sorted)
	np.savetxt("lambda_eig",eigval_sorted)
	np.save("mu_eig.npy",mu)
	np.savetxt("mu_eig",mu)

	return Q_sorted,QINV_sorted,eigval_sorted,mu

def mode_mad(traj,Q,QINV):
	import numpy as np
	import matplotlib.pyplot as plt
	import physt

	pi = np.pi
	#Get basic information from the protname.txt file
	txt = np.genfromtxt('protname.txt',dtype=str)
	protname = txt[0]
	N = int(txt[1])
	nfrs = int(txt[2])
	natoms = int(txt[3])
	nmol = int(np.loadtxt('nmol.dat'))

	#Allocate arrays containing the number of alpha-carbons per molecule
	nca = np.zeros(nmol)
	for i in range(nmol):
		nca[i] = int(np.loadtxt('nres'+str(i+1)+".dat"))
	nca_cumsum = np.cumsum(nca)

	print(protname,N,nfrs,natoms)

	traj = np.load(str(traj))
	rx = np.zeros((N,nfrs))
	ry = np.zeros((N,nfrs))
	rz = np.zeros((N,nfrs))
	lx = np.zeros((N - nmol, nfrs))
	ly = np.zeros((N - nmol, nfrs))
	lz = np.zeros((N - nmol, nfrs))
	Rinv = np.zeros((N,N))
	for numba,k in enumerate(range(0,N*nfrs,N)):
		rx[:,numba] = traj[k:k+N,0]
		ry[:,numba] = traj[k:k+N,1]
		rz[:,numba] = traj[k:k+N,2]

	#Define bond vectors
	counter = 0
	mol_counter = 0
	i = 1
	while i <= N:
		if i == nca_cumsum[mol_counter]:
			i += 1
			mol_counter += 1
			if mol_counter == nmol: break
		else:
			lx[counter,:] = rx[i,:] - rx[i-1,:] 
			ly[counter,:] = ry[i,:] - ry[i-1,:]
			lz[counter,:] = rz[i,:] - rz[i-1,:]
			counter += 1
			i += 1

	xix = np.matmul(QINV,lx)
	xiy = np.matmul(QINV,ly)
	xiz = np.matmul(QINV,lz)
	xim = np.sqrt(xix**2 + xiy**2 + xiz**2)
	xi = xix + xiy +xiz

	theta = np.arccos(xiz/xim)
	phi = np.arctan(xiy/xix)

	for a in range(N - nmol):
		for k in range(nfrs):
			if xix[a,k] <= 0.0:
				phi[a,k] += pi
			if phi[a,k] <= 0.0:
				phi[a,k] += 2*pi

	#theta = np.rad2deg(theta)
	#phi = np.rad2deg(phi)

	#Make histogram
	kT = 0.00198*float(np.loadtxt('temp'))
	fmadlist = []
	for a in range(N - nmol):
		x=xim[a,:]*np.sin(theta[a,:])*np.cos(phi[a,:])
		y=xim[a,:]*np.sin(theta[a,:])*np.sin(phi[a,:])
		z=xim[a,:]*np.cos(theta[a,:])
		h=physt.special.spherical_histogram(np.column_stack([x,y,z]),theta_bins=50,phi_bins=50,radial_bins=1)
		his=(h.densities[0]/h.densities[0].sum()).T
		fes = -kT*np.log(his)
		#fes -= fes.min()
		femax = -kT*np.log(1/nfrs)
		for j in range(fes.shape[0]):
			for k in range(fes.shape[1]):
				if np.isinf(fes[j,k]) == True:
					fes[j,k] = femax
					
		#xx,yy=np.meshgrid(np.linspace(0,180,his.shape[1]),np.linspace(0,360,his.shape[0]))
		#im=plt.contourf(xx,yy,fes,levels=np.linspace(fes.min(),fes.max(),25),cmap='gnuplot')
		#cbar=plt.colorbar(im)
		#cbar.set_label(r'Free Energy $(k_BT)$')
		#plt.contour(xx,yy,fes,levels=np.linspace(fes.min(),fes.max(),25),colors='k')
		#plt.xlabel(r'$\theta$ (deg)')
		#plt.ylabel(r'$\phi$ (deg)')
		#plt.savefig('fes'+str(a+1)+'.eps',dpi=300)
		#plt.show()
		#plt.close()
		dummy = list(np.ravel(fes))
		#Cut-off for outliers 
		cut = -kT*np.log(2/nfrs)
		for num,i in enumerate(dummy):
			if i >= cut: dummy.remove(i)
		fmad = np.median(abs(np.array(dummy) - fes.min()))

		fmadlist.append(fmad)
		np.save('fes'+str(a+1)+'.npy',fes)
		np.save('theta_phi_'+str(a+1)+'.npy',np.column_stack([xim[a,:],np.rad2deg(theta[a,:]),np.rad2deg(phi[a,:])]))
		np.savetxt("xi_"+str(a+1)+'.xvg',xi[a,:])
		np.save("xi_"+str(a+1)+'.npy',xi[a,:])
	np.savetxt('barriers_kcal.dat',np.array(fmadlist))
	np.savetxt('barriers.dat',np.array(fmadlist)/kT)

def tau_convert(eigvals,sigma,bar):
	import numpy as np
	
	#Convert kT to units of kcal/mol
	kT = 0.00198*float(np.loadtxt('temp'))
	tau = (eigvals*sigma)**-1
	tau_scaled = tau*np.exp(bar/kT)
	np.savetxt('tau.dat',np.column_stack([np.arange(1,len(tau)+1),tau]))
	np.savetxt('tau_scaled.dat',np.column_stack([np.arange(1,len(tau)+1),tau_scaled]))

def LML(Qmatrix,avbl,mu):
	import numpy as np
	
	LML = np.zeros_like(Qmatrix)
	for a in range(Qmatrix.shape[0]):
		for i in range(Qmatrix.shape[1]):
			LML[i,a] = (Qmatrix[i,a]**2*avbl**2)/mu[a]
	LML = np.sqrt(LML)

	np.savetxt("LML.dat",LML)
	np.save("LML.npy",LML)

#Need to run m1_f2py.sh first to generate the m1int 'module' file
def calc_m1(N):
	import numpy as np
	import m1int

	m1int.p2m1(N)

def calc_Rete(traj):
	import nummpy as np
	import subprocess

	traj = np.load(str(traj))
	txt = np.genfromtxt('protname.txt',dtype=str)
	protname = txt[0]
	N = int(txt[1])
	nfrs = int(txt[2])
	natoms = int(txt[3])

	print(protname,N,nfrs,natoms)

	rx = np.zeros((N,nfrs))
	ry = np.zeros((N,nfrs))
	rz = np.zeros((N,nfrs))
	Rinv = np.zeros((N,N))
	#traj = np.load(str(traj))
	for numba,k in enumerate(range(0,N*nfrs,N)):
		rx[:,numba] = traj[k:k+N,0]
		ry[:,numba] = traj[k:k+N,1]
		rz[:,numba] = traj[k:k+N,2]

	Rete = np.zeros(nfrs)
	Rete_sq = np.zeros(nfrs)

	Rete = np.sqrt((rx[-1,:] - rx[0,:])**2 + (ry[-1,:] - ry[0,:])**2 + (rz[-1,:] - rz[0,:])**2)
	Rete_sq = ((rx[-1,:] - rx[0,:])**2 + (ry[-1,:] - ry[0,:])**2 + (rz[-1,:] - rz[0,:])**2)

	np.savetxt('Rete.xvg',np.column_stack([0.2*np.arange(0,nfrs),Rete]))
	np.savetxt('Rete_sq.xvg',np.column_stack([0.2*np.arange(0,nfrs),Rete_sq]))

	subprocess.call('gmx analyze -f Rete.xvg -ac Rete_ac.xvg',shell=True)
	subprocess.call('gmx analyze -f Rete_sq.xvg -ac Retesq_ac.xvg',shell=True)

	return Rete,Rete_sq
