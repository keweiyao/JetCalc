#!/usr/bin/env python3

import numpy as np

cuts = {
'CMS':{
	'Raa':{
		'pTbins': np.array([[2,3], [3,4], [4,5], [5,6], [6,8], [8,10],
						[10,12.5], [12.5,15], [15,20], [20,25], 
						[25,30],[30,40],[40,60],[60,100]]),
		'ybins': np.array([-1.0, 1.0]),
		'cenbins':np.array([0,10],[0,100])
		  },
	'vn_HF': {
		'pTbins': np.array([[1,2],[2,3],[3,4],[4,5],[5,6],
						[6,8],[8,10],[10,15],[15,20],[20,40]]),
		'ybins': np.array([-1.0, 1.0]),
		'cenbins':np.array([0,10],[10,30],[30,50])
		 },
	'vn_ref':{
		'pTbins': np.array([0.3, 3.0]),
		'ybins': np.array([-2.4, 2.4]),
		'cenbins':np.array([0,10],[10,30],[30,50])
		}
	  },
'ALICE':{
	'Raa':{
		'pTbins': np.array([[1,2],[2,3],[3,4],[4,5],[5,6],[6,7],[7,8],[8,10],[10,12],[12,16],[16,24],[24,36], [36,50]]),
		'ybins': np.array([-0.5, 0.5]),
		'cenbins':np.array([0,10],[30,50],[60,80])
		  },
	'vn_HF': {
		'pTbins': np.array([[1,2],[2,3],[3,4],[4,5],[5,6],
							[6,7],[7,8],[8,10],[10,12],[12,16],[16,24]]),
		'ybins': np.array([-0.8, 0.8]),
		'cenbins':np.array([30,50])
		 },
	'vn_ref':{
		'pTbins': np.array([0.15, 5.0]),
		'ybins': np.array([-0.8, 0.8]),
		'cenbins':np.array([30,50])
		}
	  }
}
