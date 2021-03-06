#!/usr/bin/env python3
import numpy as np
import h5py

def fourvec_to_curvelinear(px, py, pz, E):
	pT = np.sqrt(px**2+py**2)
	phi = np.arctan2(py, px)
	y = .5*np.log((E+pz)/(E-pz))
	M = np.sqrt(E**2 - px**2 - py**2 - pz**2)
	return M, pT, y, phi

def threevec_to_curvelinear(px, py, pz, M):
	pT = np.sqrt(px**2+py**2)
	phi = np.arctan2(py, px)
	E = np.sqrt(M**2 + px**2 + py**2 + pz**2)
	y = .5*np.log((E+pz)/(E-pz))
	return E, pT, y, phi

def EventPlane(particle_dict, pTcut, ycut, order=6):
	pT = particle_dict['pT']
	phi = particle_dict['phi']
	y = particle_dict['y']

	results = {'Qn':{}, 'M':{}}
	cut = (pTcut[0]<pT) & (pT<pTcut[1]) & (ycut[0]<y) & (y<ycut[1])
	phi = phi[cut]
	for n in range(1, 1+order):
		results['Qn'][n] = np.exp(1j*n*phi).sum()
		results['M'][m] = cut.sum()
	return results

def Qvector(particle_dict, pTbins, ybins, pid_POI, order=4):
	pT = particle_dict['pT']
	phi = particle_dict['phi']
	y = particle_dict['y']
	if particle_dict['w'] is None:
		w = np.ones_like(pT)
	else:
		w = particle_dict['w']

	# calculate Qn
	results = {pid : 
				{'M': np.zeros([len(pTbins), len(ybins)], dtype=np.float),
				 'Qn': np.zeros([len(pTbins), len(ybins), order], dtype=np.complex)}
				for pid in pid_POI}
	for pid in pid_POI:
		pid_cut = (pid == particle_dict['pid'])
		for i, (pl, ph) in enumerate(pTbins):
			pT_cut = (pl<=pT) & (pT<=ph)
			for j, (yl, yh) in enumerate(ybins):
				y_cut = (yl<=y) & (y<=yh)			
				indices = pid_cut & pT_cut & y_cut
				# apply cut
				phi_POI = phi[indices]
				w_POI = w[indices]
		
				for n in range(1, 1+order):
					results[pid]['Qn'][i,j,n-1] = np.sum(np.exp(1j*n*phi_POI)*w_POI)
				results[pid]['M'][i,j] = w_POI.sum()
	return results

def Yield(particle_dict, pTbins, ybins, pid_POI, order=4):
	pT = particle_dict['pT']
	phi = particle_dict['phi']
	y = particle_dict['y']
	if particle_dict['w'] is None:
		w = np.ones_like(pT)
	else:
		w = particle_dict['w']

	# calculate Qn
	results = {pid : np.zeros([len(pTbins), len(ybins)], dtype=np.float) 
					for pid in pid_POI}
	for pid in pid_POI:
		pid_cut = (pid == particle_dict['pid'])
		for i, (pl, ph) in enumerate(pTbins):
			pT_cut = (pl<=pT) & (pT<=ph)
			for j, (yl, yh) in enumerate(ybins):
				y_cut = (yl<=y) & (y<=yh)			
				indices = pid_cut & pT_cut & y_cut
				# apply cut
				w_POI = w[indices]
				results[pid][i,j] = w_POI.sum()
	return results

