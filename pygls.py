#!/usr/bin/env python

import os
import sys
import getopt
import time
import nibabel as nib
import numpy as np
import math
import shutil

############################## FUNCTIONS ######################################
def transformToFreqDomain(tlen,autocorr_len,func_data,acp,design_mat):
	"""Transforms temporal data into the frequency domain.

	All temporal signals are zero-padded to length:
		tlen + 2*autocorr_len + 2
	This avoids any wrap-around artifacts that would occur from convolving the
	auto-correlation filter with the functional data.

	Args:
		tlen: number of TR's in the functional data. When func_data is present,
			tlen = func_data.shape[0]
		autocorr_len: maximum lag of the auto-correlation filter. When acp is 
			present, autocorr_len = acp.shape[0]
		func_data: 2D numpy array of functional data.  Each column is a 
			functional time-series associated with a given voxel.  A 1D FFT
			is computed along the rows (time dimension).
			If empty, this input is ignored.
		acp: 2D numpy array of auto-correlation estimates. Each column contains
			noise auto-correlation estimates for a given voxel.  A 1D FFT 
			is computed along the rows (time-lag dimension).
			If empty, this input is ignored.
		design_mat: 2D numpy array (the design matrix). Each column corresponds
			to a regressor in the model.  A 1D FFT is computed along the rows
			(time dimension).
			If empty, this input is ignored.

	Returns:
		A 3-element tuple containing the FFT coefficients for the func_data, 
		acp, and design_mat. If an input is empty, the corresponding FFT output
		will also be empty.  
	"""

	eps = 1e-5	
	FFT_N = int(tlen+2*autocorr_len+2)
	if func_data.size > 0:
		if len(func_data.shape) == 2:
			func_data_mean = np.mean(func_data,0)
			FUNC_DATA = np.fft.fft(func_data - func_data_mean, n=FFT_N, axis=0)
		else:
			FUNC_DATA = []
	else:
		FUNC_DATA = []

	if acp.size> 0:
		acp_sz = acp.shape
		if len(acp_sz) == 2:
			corr_filt = np.zeros([FFT_N,acp_sz[1]])
			corr_filt[0:autocorr_len,:] = acp
			corr_filt[::-1,:][0:autocorr_len-1] = corr_filt[1:autocorr_len,:] 
			corr_filt = corr_filt - np.mean(corr_filt,0)
			ACP = np.real(np.fft.fft(corr_filt,n=FFT_N,axis=0))

			locs = np.argwhere(ACP>eps)
			locs = np.transpose(locs)
			ACP[locs[0],locs[1]] = 1/np.sqrt(ACP[locs[0],locs[1]])
			ACP = ACP/np.tile(np.sqrt(sum(ACP**2,0)/FFT_N),[FFT_N,1])
		else:
			ACP = []
	else:
		ACP = []

	if design_mat.size > 0:
		design_mat_mean = np.mean(design_mat,0)
		DESIGN_MAT = np.transpose(np.fft.fft(design_mat - design_mat_mean, 
											 n=FFT_N, 
											 axis=0))
	else:
		DESIGN_MAT = []
	return (FUNC_DATA,ACP,DESIGN_MAT)
	

def setup_susan_filter(func_data,vox_msk_inds,full_example_slice,
					   susan_ms,dims,shape):
	"""Precomputes information used by spatial SUSAN filtering.

	Args:
		func_data:
		vox_msk_inds:
		full_example_slice:
		susan_ms:
		dims:
		shape:

	Returns:
		A dictionary containing pre-computed information for the SUSAN filter
		Example:
			{'mask_sizes':,
			 'vox_sub2ind':,
			 'vox_coords':,
			 'vox_in_mask_inds':,
			 'I':,
			 'I2':
			}
	"""

	sx,sy,sz = shape[0:3]
	xdim,ydim,zdim = dims[0:3]

	N = len(full_example_slice)
	msk_sz = len(vox_msk_inds)

	example_slice = full_example_slice[vox_msk_inds]

	# Now compute susan_mode and susan_thresh based on example_slice
	x1,x2 = np.histogram(example_slice,max(1,math.floor(msk_sz/200)))
	x1max = max(x1)
	ind = np.nonzero(x1==x1max)
	print ind
	if len(ind[0]) > 1:
		susan_mode = x2[ind[0][0]+1]
	else:
		susan_mode = x2[ind[0]+1]
	#print 'mode = %f' % susan_mode
	example_slice[np.where(example_slice < susan_mode)] = np.float32(susan_mode)
	susan_sig = math.floor(math.sqrt(sum((example_slice-susan_mode)**2)/msk_sz))
	#print 'sig = %d' % int(susan_sig)

	## Now setup the susan structure for spatial filtering
	susan_thresh = math.floor(susan_sig/3)
	susan = dict()
	susan['mask_sizes'] = 2*np.ceil(2*susan_ms/np.array([xdim,ydim,zdim]))+1
	susan['vox_sub2ind'] = np.reshape(
							np.arange(0,N,dtype=np.int),[sx,sy,sz],order='F')
	susan['vox_coords'] = np.array(
							np.unravel_index(
								np.arange(0,N),[sx,sy,sz],order='F'))
	susan['vox_in_mask_inds'] = -1*np.ones(N,dtype=np.int)
	susan['vox_in_mask_inds'][vox_msk_inds] = np.arange(0,msk_sz)
	# Now scale the dimensions so that the weights will be just
	# exp(-||susan.I(:,i)-susan.I(:,j)||^2)
	Iscale = np.hstack((1.0/susan_thresh, 
						np.array([xdim,ydim,zdim])/(math.sqrt(2)*susan_ms)))
	# Scale the rows of susan['I']
	susan['I'] = Iscale[:,np.newaxis]*np.vstack((full_example_slice, 
												 susan['vox_coords']))
	# Store sum of squares for faster distance computation later
	susan['I2'] = np.sum(susan['I']**2,axis=0)

	return susan


def susan_filt(data,s):
	"""Perform SUSAN spatial filtering of a given set of data.  Implementation
	matches FSL's version:
	http://fsl.fmrib.ox.ac.uk/fsl/fslwiki/SUSAN

	Args:
		data: a 2D numpy array.  Each column is a time-series associated with a 
			given voxel.
		s: dictionary containing information about the smoothing filter. 

	Returns:
		data_smoothed:
	"""

	sx,sy,sz = s['vox_sub2ind'].shape
		
	vox_sub2ind = s['vox_sub2ind']
	vox_coords = s['vox_coords']
	I = s['I']
	I2 = s['I2']
		  
	msx,msy,msz = (s['mask_sizes']-1)/2
	
	data_smoothed = np.zeros(data.shape,order='F')
	
	for vox_ind in xrange(I.shape[1]):
		v2 = s['vox_in_mask_inds'][vox_ind]
		if v2 < 0:
			continue

		vc = vox_coords[:,vox_ind]
		# Determine the neighboring voxels
		nbr_vox_inds = np.reshape(
						vox_sub2ind[max(vc[0]-msx,0):min(vc[0]+msx,sx)+1,
									max(vc[1]-msy,0):min(vc[1]+msy,sy)+1,
									max(vc[2]-msz,0):min(vc[2]+msz,sz)+1],-1)
		
		#Now compute sum-of-squared distances between voxel and its neighbors
		wts = np.exp(2*(I[:,nbr_vox_inds].T.dot(I[:,vox_ind])) 
						- (I2[vox_ind] + I2[nbr_vox_inds]))

		inds = s['vox_in_mask_inds'][nbr_vox_inds]
		locs = np.where(inds >= 0)
		data_smoothed[:,v2] = (data[:,inds[locs]].dot(wts[locs]))/np.sum(wts)

	return data_smoothed

	
def computeTemporalAutoCorrEst(func_data,design_mat,susan,autocorr_len):
	"""Estimate noise auto-covariance at each voxel.

	Args:
		func_data: 2D numpy array containing the functional data. Each column
			contains the functional time-series from a given voxel.
		design_mat: 2D numpy array containing the design matrix.  Columns 
			correspond to regressors.
		susan: a dictionary of pre-computed information created by 
			setup_susan_filter.
		autocorr_len: maximum lag of noise auto-covariance to estimate

	Returns:
		corr_est: A 2D numpy array containing the auto-covariance estimates 
			at each voxel.
			corr_est[i,j] = noise auto-covariance estimate for voxel j at lag=i
			E.g., if n is an estimate of the noise at voxel j, then
			corr_est[i,j] = E[n[k]*n[k+i]]
			Here, i is in the range [0,autocorr_len-1]
	"""

	tlen = func_data.shape[0]
	# De-mean the data and design
	func_data = func_data - np.mean(func_data,axis=0)
	design_mat = design_mat - np.mean(design_mat,axis=0) 
	# Estimate regression coefficients assuming white noise
	# Compute residuals assuming iid noise
	# This will be used as an esetimate for the noise at each voxel
	resids = func_data - design_mat.dot(
							np.linalg.solve(
								design_mat.T.dot(design_mat),
								design_mat.T.dot(func_data)))
	
	# Compute PSD of residuals
	F = np.abs(np.fft.fft(resids, n=int(2*tlen+1), axis=0))**2
	corr_est = np.real(np.fft.ifft(F, axis=0)[0:autocorr_len,:])	
	del F
	# Now scale the corr_est estimates
	corr_est = (1.0/(tlen - np.arange(1,autocorr_len+1)))[:,np.newaxis]*(corr_est*(1.0/(np.sum(resids**2,axis=0)/(tlen-1))))
	# Spatially smooth the auto-corr estimates	
	# For whatever reason, FSL does not smooth the lag=0 estimates
	corr_est[1:,:] = susan_filt(corr_est[1:,:],susan)
	
	# Now window using a Tukey window
	# It should really be 0.5*(1+cos(pi*(0:autocorr_len-1)'/autocorr_len))
	# But it is set to the following below to match FSL
	tukey_win = 0.5*(1+np.cos(math.pi*np.arange(1,autocorr_len+1)/autocorr_len))
	# Window the auto-corr
	corr_est = tukey_win[:,np.newaxis]*corr_est
	return corr_est


def readDesignFSFFile(designfsffile):
	"""Extract information from FSL design.fsf file

	Args:
		designfsffile: string -- full path to design.fsf file

	Returns: 
		A dictionary containing relevant information extracted from the FSF 
		file.  For example:
		{'designmatfile': 'design.mat',
		 'designconfile': 'design.con',
		 'func_datafile': 'filtered_func_data.nii.gz',
		 'mask_datafile': 'mask.nii.gz',
		 'outputdirs': ['~/my_output_dir/'],
		 'susan_ms': 5,
		 'threshac1files': ['__ESTIMATE__']}
	"""

	# First, we must parse the designfsffile
	fsf = {'full':designfsffile}
	drive,designfsffile = os.path.splitdrive(designfsffile)
	if drive != '':
		raise Exception('This script has not been tested on Windows-like machines.')
	fsf['path'],fsf['filename'] = os.path.split(designfsffile)
	fsf['filename_noext'],fsf['extension'] = os.path.splitext(fsf['filename'])
	fsf['full_noext'] = ''.join([fsf['path'],'/',fsf['filename_noext']])

	fp = open(fsf['full'])
	fplist = fp.readlines()
	fp.close()

	# FSL hard-codes 5mm smoothing of auto-correlation estimates
	# See featlib.tcl in /usr/local/fsl/src/featlib.tcl
	susan_ms = 5
	threshac1file = ''
	outputdir = ''
	confoundev = ''
	featfile = ''
	for line in fplist:
		if line[0:18] == 'set fmri(analysis)':
			b = line.rstrip().split()
			if len(b) != 3:
				raise Exception('Could not parse fmri analysis line')
			analysis_num = int(b[2])
			if analysis_num != 2 and analysis_num != 6:
				raise Exception('Cannot run pre-stats')
		if line[0:22] == 'set fmri(prewhiten_yn)':
			b = line.rstrip().split()
			if len(b) != 3:
				raise Exception('Could not parse prewhitening option')			
			prewhiten = int(b[2])			
			if prewhiten:
				threshac1file = '__ESTIMATE__'
			else:
				threshac1file = '__NO_ESTIMATE__'				
		if line[0:16] == 'set fmri(smooth)':
			pass
			#b = line.rstrip().split()
			#if len(b) != 3:
			#	raise Exception('Could not parse spatial smoothing FWHM parameter')
			#susan_ms = int(b[2])
		if line[0:19] == 'set fmri(outputdir)':
			b = line.split('"')
			# Should return something like this:
			#	b[0] = "set feat_files(1) "
			#	b[1] = "/path/to/output_directory/"
			#	b[2] = ""
			if len(b) != 3:
				raise Exception('Could not parse output directory from design.fsf file')
			outputdir = b[1]
		elif line[0:23] == 'set confoundev_files(1)':
			b = line.split('"')
			# Should return something like this:
			#	b[0] = "set confoundev_files(1) "
			#	b[1] = "/path/to/confoundev_file/confound.par"
			#	b[2] = ""
			if len(b) != 3:
				raise Exception('Could not parse confound ev file')
			confoundev = b[1]
		elif line[0:17] == 'set feat_files(1)':
			b = line.split('"')
			# Should return something like this:
			#	b[0] = "set feat_files(1) "
			#	b[1] = "/path/to/feat_file/filtered_func_data.nii.gz"
			#	b[2] = ""
			if len(b) != 3:
				raise Exception('Could not parse feat_file')
			featfile = b[1]

	if len(threshac1file) == 0:
		raise Exception('Could not parse prewhitening option')
	if susan_ms < 0:
		raise Exception('Could not parse spatial smoothing FWHM parameter')

	# Create the output directory if not already present	
	if not os.path.exists(outputdir):
		os.makedirs(outputdir)

	# Check to see if there is a recognizable extension on the featfile
	if featfile[-6:] != 'nii.gz':
		featfile = featfile + '.nii.gz'

	# feat_model puts the design/contrast files in the location of the fsf file
	# so copy it to the output directory
	assert_design_filename = 1
	if assert_design_filename:
		fsf['filename'] = 'design.fsf'
		fsf['filename_noext'] = 'design'

	shutil.copyfile(fsf['full'],'%s/%s'%(outputdir,fsf['filename']))

	# Now we need to call feat_model
	exec_str = '%s/bin/feat_model %s/%s %s' % (os.environ['FSLDIR'], 
											   outputdir, 
											   fsf['filename_noext'], 
											   confoundev)
	ret = os.system(exec_str)
	if ret != 0:
		raise Exception('Could not run feat_model')

	# Now we need to create the mask
	exec_str = '%s/bin/fslmaths %s -Tmin -bin %s/mask_tmp -odt char' % (os.environ['FSLDIR'], featfile, outputdir)
	ret = os.system(exec_str)
	if ret != 0:
		raise Exception('Could not run fslmaths')

	# Extract the minimum intensity relative to the mask
	exec_str = '%s/bin/fslstats %s -k %s/mask_tmp -R | /usr/bin/awk \'{ print }\' -' % (os.environ['FSLDIR'],featfile,outputdir)
	ret = os.popen(exec_str).read()
	binarise_thresh_min = float(ret.split()[0])

	# Now we can delete the mask_tmp
	os.remove('%s/mask_tmp.nii.gz' % outputdir)

	# Now re-create the mask based on mean and variance of the filtered-func image
	exec_str = '%s/bin/fslmaths %s -Tmean -thr %f -bin %s/tmp1' % (os.environ['FSLDIR'],featfile,binarise_thresh_min,outputdir)
	ret = os.system(exec_str)
	if ret != 0:
		raise Exception('Could not run fslmaths')
	
	exec_str = '%s/bin/fslmaths %s -Tstd -thr 0.00001 -bin -mul %s/tmp1 %s/mask' % (os.environ['FSLDIR'],featfile,outputdir,outputdir)
	ret = os.system(exec_str)
	if ret != 0:
		raise Exception('Could not run fslmaths')

	# Now we can delete the temporary file
	os.remove('%s/tmp1.nii.gz' % outputdir)

	return {'designmatfile' : '%s.mat' % fsf['filename_noext'],
			'designconfile' : '%s.con' % fsf['filename_noext'],
			'func_datafile' : featfile,
			'mask_datafile' : '%s/mask.nii.gz' % (outputdir),
			'outputdirs' : [outputdir],
			'susan_ms' : susan_ms,
			'threshac1files' : [threshac1file]
			}


def generateDesignMatrixFromTiming(stimTimes,TR,hrf='Double-Gamma'):
	"""Generates design matrix from stimulus timing information

	Args:
		stimTimes: (units of seconds)
		TR: (units of seconds)
		hrf:

	Returns:
		design_mat: 
	"""

	pass


def readDesignAndContrastFiles(designmatfile,designconfile):
	"""Read design matrix and contrasts from file.

	Args:
		designmatfile: string -- filename of design.mat file.  This file is 
			created by FSL's feat_model and stores the design matrix
		designconfile: string -- filename of design.con file.  This file is 
			created by FSL's feat_model and stores the contrasts between 
			regressors.

	Returns:
		A two-element tuple containing the design matrix and contrasts.
		For example:
			(design_mat, contrasts)
		design_mat is a 2D numpy array. Each column is a regressor
		contrasts is a 2D numpy array.  Each column specifies a contrast between
		regressors.  For example, if the first contrast takes the difference
		between the first and second regressors, then contrasts[:,1] would be:
			contrasts[:,1] = np.array([1,-1,0,0,...,0])
	"""

	fp = open(designmatfile)
	fplist = fp.readlines()
	fp.close()
	C = fplist.pop(0)
	numWaves = int(C.split('\t')[1].rstrip())
	C = fplist.pop(0)
	numPoints = int(C.split('\t')[1].rstrip())#replace('\n',''))
	while fplist.pop(0)[0:7] != '/Matrix':
		pass

	C = fplist
	if len(C) != numPoints:
		raise Exception("Could not parse design.mat file")

	design_mat = np.zeros([numPoints,numWaves])
	for i in xrange(len(C)):
		design_mat[i,:] = np.asarray([float(_) for _ in C[i].rstrip().split()])
		
	nmp = design_mat.shape[1]

	fp = open(designconfile)
	fplist = fp.readlines()
	fp.close()
	for i in xrange(len(fplist)):
		C = fplist.pop(0)
		if C[:13] == '/NumContrasts':
			numContrasts = int(C.split('\t')[1].rstrip())
			break

	while fplist.pop(0)[0:7] != '/Matrix':
		pass

	C = fplist
	if len(C) != numContrasts:
		raise Exception("Could not parse design.con file")

	contrasts = np.zeros([nmp,numContrasts])
	for i in xrange(len(C)):
		contrasts[:,i] = np.asarray([float(_) for _ in C[i].rstrip().split()])

	return(design_mat,contrasts)


def getDefaultSaveOpts():
	return {'cope':1,
			'dof':1,
			'neff':0,
			'pe':0,
			'threshac1':0,
			'tstat':0,
			'varcope':1}


def run_pygls(fsffile,saveopts=dict()):
	"""

	Args:
		fsffile:
		saveopts:

	Returns:
		Nothing. All results are saved in the output directory specified in the
		design.fsf file.  File output format matches that of FSL.

	"""

	tstart = time.time()

	defSaveOpts = getDefaultSaveOpts()
	for key in defSaveOpts.keys():
		if key not in saveopts:
			saveopts[key] = defSaveOpts[key]

	inpmode = 'fsf'
	if inpmode == 'fsf':
		ret = readDesignFSFFile(fsffile)
		func_datafile = ret['func_datafile']
		mask_datafile = ret['mask_datafile']
		outputdirs = ret['outputdirs']
		threshac1files = ret['threshac1files']
		designmatfile = ret['designmatfile']
		designconfile = ret['designconfile']

		susan_ms = ret['susan_ms']
	else:
		raise Exception('Unknown input mode')

	## Read in filtered func data
	func_nii = nib.load(func_datafile)
	func_data = func_nii.get_data()
	scales = func_nii.get_affine()
	shape = func_nii.shape
	sx,sy,sz,tlen = shape[0:4]
	del func_nii
	dims = np.around(np.abs(scales.diagonal(0)))
	xdim,ydim,zdim = dims[0:3]

	## Read in the mask file
	mask_nii = nib.load(mask_datafile)
	ms = mask_nii.shape
	mask_data = mask_nii.get_data()
	vox_msk_inds = np.nonzero(mask_data)
	vox_msk_inds = np.array(vox_msk_inds)
	msk_sz = np.shape(vox_msk_inds)[1]
	vox_msk_inds = np.sort(
					np.ravel_multi_index(
						[vox_msk_inds[2],vox_msk_inds[1],vox_msk_inds[0]],
						dims=(ms[2],ms[1],ms[0])))
	vrm = np.unravel_index(vox_msk_inds,[sx,sy,sz],order='F')	
	del mask_nii, mask_data

	autocorr_len = int(math.floor(math.sqrt(tlen)))
	
	## Now re-shape functional data based on voxel mask
	func_data = np.reshape(np.transpose(func_data,[3,0,1,2]),
						   (tlen,sx*sy*sz),
						   order='F')
	N = func_data.shape[1]
	full_example_slice = func_data[math.floor(tlen/2)-1,:]
	func_data = func_data[:,vox_msk_inds]

	susan = setup_susan_filter(func_data,
							   vox_msk_inds,
							   full_example_slice,
							   susan_ms,
							   dims,
							   shape)
	del full_example_slice

	## Transform the functional data into the frequency domain
	[FUNC_DATA,_,_] = transformToFreqDomain(tlen,
											autocorr_len,
											func_data,
											np.array([]),
											np.array([]))
	FFT_N = FUNC_DATA.shape[0]
	sqrtFFT_N = math.sqrt(FFT_N)

	# Frequency-domain method:
	# We need to compute a portion of the DFT matrix to account for the zero padding
	twid = complex(0,1)*(2*math.pi/FFT_N)
	DFTV = np.exp(twid * np.outer(
							np.arange(tlen,FFT_N),np.arange(0,FFT_N)))/sqrtFFT_N
	DFTWide = np.hstack((np.real(DFTV), -1*np.imag(DFTV)))
	del DFTV, twid

	FTall_orig = np.vstack((np.real(FUNC_DATA), np.imag(FUNC_DATA)))/sqrtFFT_N
	del FUNC_DATA
	for dNum in xrange(len(outputdirs)):
		outputdir = outputdirs[dNum]    
		threshac1file = threshac1files[dNum]
		outstatsdir = '%s/stats/' % outputdir
		if os.path.exists(outstatsdir):
			shutil.rmtree(outstatsdir)

		os.makedirs(outstatsdir)
		designmatfilecurr = '%s/%s' % (outputdir, designmatfile)
		designconfilecurr = '%s/%s' % (outputdir, designconfile)

		design_mat,contrasts = readDesignAndContrastFiles(designmatfilecurr,
											designconfilecurr)
		del designmatfilecurr, designconfilecurr
		nmp = design_mat.shape[1]
		numContrasts = contrasts.shape[1]

		hasAutoCorrEst = 1
		# Read in the auto-corr file
		if threshac1file == '__ESTIMATE__':
			# Do the temporal auto-correlation here
			acp = computeTemporalAutoCorrEst(func_data,
											 design_mat,
											 susan,
											 autocorr_len)
			if dNum == len(outputdirs) - 1:
				del func_data, susan
			
			if saveopts['threshac1'] == 1:
				# Save out the acp to file
				ACP_SV = np.zeros([sx,sy,sz,autocorr_len])
				ACP_TMP = np.zeros([sx,sy,sz])
				for i in xrange(autocorr_len):
					ACP_TMP[vrm] = acp[i,:]
					ACP_SV[:,:,:,i] = ACP_TMP
				pe_img = nib.Nifti1Image(ACP_SV,scales)
				nib.save(pe_img, '%s/threshac1.nii.gz' % (outstatsdir))
				del ACP_SV, ACP_TMP, pe_img
		elif threshac1file == '__NO_ESTIMATE__':
			# No temporal autocorrelation filtering
			hasAutoCorrEst = 0
			acp = 0
		else:
			acp = nib.load(thresac1file)
			if acp[3].size != autocorr_len:
				print 'Problem here!'
				acp = np.reshape(np.transpose(acp,[3,0,1,2]),
								 (autocorr_len,sx*sy*sz))
				acp = acp[:,vox_msk_inds]

		if hasAutoCorrEst == 0:
			# There is a much faster way of doing all this
			print('There is a faster way of doing this!')
			# Do some stuff here and then
			ACP = np.ones([FFT_N,msk_sz])
			# continue
			[_,_,DESIGN_MAT] = transformToFreqDomain(tlen,
													 autocorr_len,
													 np.array([]),
													 np.array([]),
													 design_mat)
		else:
			# Now transform the auto-corr estimates into the frequency domain
			[_,ACP,DESIGN_MAT] = transformToFreqDomain(tlen,
													   autocorr_len,
													   np.array([]),
													   acp,
													   design_mat)

		ACP2 = np.vstack((ACP,ACP))
		DMTall = np.transpose(
					np.hstack(
					   (np.real(DESIGN_MAT),np.imag(DESIGN_MAT))))/sqrtFFT_N
		del ACP, acp, design_mat, DESIGN_MAT

		# Now filter the functional data
		FTall = ACP2*FTall_orig
		FUNC_DATA_CORRECTION = DFTWide.dot(FTall)
		ftf = (np.transpose(np.sum(FTall**2,axis=0) 
					- np.sum(FUNC_DATA_CORRECTION**2, axis=0)))
		# Compute dof and save t file
		dof = tlen - nmp
		if saveopts['dof'] == 1:
			fp = open('%s/dof' % outstatsdir,'w+')
			fp.write('%d\n' % dof)
			fp.close()

		betas_msk = np.zeros([nmp,msk_sz])
		neffinv = np.zeros([numContrasts,msk_sz])
		sigmasquared_msk = np.zeros([msk_sz])
		ctmp = np.hstack((np.zeros([nmp,1]), contrasts))
		for j in xrange(msk_sz):
			DM = (DMTall.T*ACP2[:,j]).T
			DESIGN_MAT_CORRECTION = DFTWide.dot(DM)

			mtm = DM.T.dot(DM) - DESIGN_MAT_CORRECTION.T.dot(DESIGN_MAT_CORRECTION)
			mtf = DM.T.dot(FTall[:,j]) - DESIGN_MAT_CORRECTION.T.dot(FUNC_DATA_CORRECTION[:,j])
			ctmp[:,0] = mtf
			c = np.linalg.solve(mtm,ctmp)

			betas_msk[:,j] = c[:,0]
			neffinv[:,j] = np.sum(contrasts*c[:,1:],axis=0)

			sigmasquared_msk[j] = (ftf[j] - c[:,0].dot(2*mtf - mtm.dot(c[:,0])))/dof
					
		c = np.zeros([sx,sy,sz])
		if saveopts['pe'] == 1:
			for j in xrange(nmp):
				filenum = j+1
				c[vrm] = betas_msk[j,:]
				pe_img = nib.Nifti1Image(c, scales)
				nib.save(pe_img, '%s/pe%i.nii.gz' % (outstatsdir,j+1))

		if saveopts['cope'] == 1:
			copes_msk = contrasts.T.dot(betas_msk)
		del betas_msk
		if saveopts['varcope'] == 1:
			varcopes_msk = neffinv*sigmasquared_msk
		if saveopts['neff'] == 1:
			neff_msk = 1.0/neffinv
		del neffinv

		for j in xrange(numContrasts):
			if saveopts['neff'] == 1:
				c[vrm] = neff_msk[j,:]
				pe_img = nib.Nifti1Image(c, scales)
				nib.save(pe_img, '%s/neff%d.nii.gz' % (outstatsdir,j+1))
			
			if saveopts['cope'] == 1:
				c[vrm] = copes_msk[j,:]
				pe_img = nib.Nifti1Image(c, scales)
				nib.save(pe_img, '%s/cope%d.nii.gz' % (outstatsdir,j+1))
			
			if saveopts['varcope'] == 1:
				c[vrm] = varcopes_msk[j,:]
				pe_img = nib.Nifti1Image(c, scales)
				nib.save(pe_img, '%s/varcope%d.nii.gz' % (outstatsdir,j+1))
			
			if saveopts['tstat'] == 1:
				vctmp = np.sqrt(varcopes_msk[j,:])
				vctmp[np.where(vctmp<1e-5)] = 1				
				c[vrm] = copes_msk[j,:]/vctmp
				pe_img = nib.Nifti1Image(c, scales)
				nib.save(pe_img,'%s/tstat%d.nii.gz' % (outstatsdir,j+1))
		
	elapsed_time = time.time() - tstart
	print "pygls completed...took %.2f seconds" % elapsed_time


if __name__ == '__main__':
	defSaveOpts = getDefaultSaveOpts()

	usage_str = """pygls.py -f <design.fsf> [options]
	Available options:
		--cope=[0|1]:        Save contrast parameter estimates (copeXX.nii.gz) (DEFAULT: %d)
		--dof=[0|1]:         Save degrees of freedom file (dof) (DEFAULT: %d)
		--fsf=<design.fsf>:  Equivalent to -f <design.fsf>
		-h:                  Display usage and exit
		--neff=[0|1]:        Save neff images (neffXX.nii.gz) (DEFAULT: %d)
		--pe=[0|1]:          Save parameter estimates (peXX.nii.gz) (DEFAULT: %d)
		--threshac1=[0|1]:   Save auto-correlation estimates (threshac1.nii.gz) (DEFAULT: %d)
		--tstat=[0|1]:       Save tstat images (tstatXX.nii.gz) (DEFAULT: %d)
		--varcope=[0|1]:     Save covariance of copes (varcopeXX.nii.gz) (DEFAULT: %d)""" % (defSaveOpts['cope'],
		   defSaveOpts['dof'],
		   defSaveOpts['neff'],
		   defSaveOpts['pe'],
		   defSaveOpts['threshac1'],
		   defSaveOpts['tstat'],
		   defSaveOpts['varcope'])

	saveopts = {}
	argv = sys.argv[1:]
	# Parse the command line arguments
	try:
		opts, args = getopt.getopt(argv,'hf:',['cope=',
											   'dof=',
											   'fsf=',
											   'neff=',
											   'pe=',
											   'tstat=',
											   'varcope='])
	except getopt.GetoptError:
		print getopt
		print usage_str,
		sys.exit()

	if not len(opts):
		print usage_str,
		sys.exit()
	for opt,arg in opts:
		if opt == '-h':
			print usage_str,
			sys.exit()
		elif opt in ('-f','--fsf'):
			inpmode = 'fsf'
			fsffile = arg
		elif opt in ('--pe','--cope','--varcope','--dof','--neff','--tstat'):
			saveopts[opt[2:]] = arg
		
	del usage_str, defSaveOpts
	run_pygls(fsffile,saveopts)
