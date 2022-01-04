# SSPACE-MS analysis
#
# Authors:          Carmen Tur Gomez, University College London
#                  <c.tur@ucl.ac.uk> <carmen.tur.gomez@gmail.com>
#                   Francesco Grussu, University College London
#		    <f.grussu@ucl.ac.uk> <francegrussu@gmail.com>
#
# Code released under BSD Two-Clause license
#
# Copyright (c) 2021-2022 University College London. 
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
# 
# 1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
# 2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
# 
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
# 
# The views and conclusions contained in the software and documentation are those
# of the authors and should not be interpreted as representing official policies,
# either expressed or implied, of the FreeBSD Project.


### Load useful modules
import argparse
import nibabel as nib
import numpy as np
import sys
import warnings
from skimage import measure


# Function to calculate SSPACE-MS metrics
def getSSPACEMSmetrics(mask_data,mask_resmm,refvox):
	''' Get SSPACE-MS metrics on a binary lesion mask provided as a numpy array

	    INTERFACE
	    metrics = getSSPACEMSmetrics(mask_data,mask_resmm,refvox)

	    PARAMETERS
		- mask_data:  3D numpy array storing the lesion mask (1 in lesional voxel, 0 in non-lesional voxels)
		- mask_resmm: 3-element numpy array storing the resolution in mm along the 3 image dimensions (i, j, k)
		- refvox:     3-element numpy array storing the position in image coordinate of the reference voxel 
                             for centre-of-mass calculations (i.e. position as (i, j, k), example: refvox=[0,2,10])  

	    OUTPUT
		- metrics: array of 17 elements storing SSPACE-MS metrics obtained from 
			   S = cov(voxel_position) (covariance matrix of voxel positions) 
			   and COM(mask_data) as

				metrics[0] --> axial_cov: L1 or first eigenvalue (axial covariance) [mm^2]
				metrics[1] --> second_eig: L2 or second eigenvalue [mm^2] 
				metrics[2] --> third_eig: L3 or third eigenvalue [mm^2]
				metrics[3] --> radial_cov: radial covariance [mm^2]
				metrics[4] --> mean_cov: mean covariance [mm^2]
				metrics[5] --> fa_cov: covfa or covariance fractional anisotropy (ranging in [0; 1], with 1 maxiumum 
							anisotropy)

						fa_cov = (sqrt(3)/sqrt(2))*
							  sqrt((axial_cov - mean_cov)^2 + ...
							  (second_eig - mean_cov)^2 + ...
							  (third_eig - mean_cov)^2 ) / sqrt( axial_cov^2+second_eig^2+third_eig^2)

				metrics[6] --> plan_cov: covpl or covariance planarity (ranging in [0; 1], with 1 maximum planarity)
							  it is defined as
								
							  plan_cov = (2/3)*(L2 - L3)/MC
							  
							  where MC = (L1 + L2 + L3)/3 is the mean covariance
								
				metrics[7] --> sph_cov: covsph or covariance sphericity (ranging in [0; 1], with 1 maximum sphericity)
							  it is defined as
								
							  sph_cov = L3/MC
							  
							  where MC = (L1 + L2 + L3)/3 is the mean covariance
				                                                         
				metrics[8] --> tot_les_load: lesld or total lesion load (mm^3)
				metrics[9] --> pos_cov[0,0]: sxx, element of the covariance matrix [mm^2] 
					       (i,j,k are image space coordinates converted in mm from spatial resolution)
				metrics[10] --> pos_cov[1,1]: syy, element of the covariance matrix [mm^2]
					       (i,j,k are image space coordinates converted in mm from spatial resolution)
				metrics[11] --> pos_cov[2,2]: szz, element of the covariance matrix [mm^2]
					       (i,j,k are image space coordinates converted in mm from spatial resolution)
				metrics[12] --> pos_cov[0,1]: sxy, element of the covariance matrix [mm^2]
					       (i,j,k are image space coordinates converted in mm from spatial resolution)
				metrics[13] --> pos_cov[0,2]: sxz, element of the covariance matrix [mm^2]
					       (i,j,k are image space coordinates converted in mm from spatial resolution)
				metrics[14] --> pos_cov[1,2]: syz , element of the covariance matrix [mm^2]
					       (i,j,k are image space coordinates converted in mm from spatial resolution)
				metrics[15] --> COMi: first coordinate of the centre of mass of the mask with resppect to the 
                                                     reference voxel [mm] (COMx = res_i*(i_COM - i_reference)
				metrics[16] --> COMj: second coordinate of the centre of mass of the mask with resppect to the 
                                                     reference voxel [mm] (COMy = res_j*(j_COM - j_reference)
				metrics[17] --> COMk: third coordinate of the centre of mass of the mask with resppect to the 
                                                     reference voxel [mm] (COM1 = res_k*(k_COM - k_reference)

	    DEPENDENCIES		
		nibabel, numpy, skimage

	    Authors: Carmen Tur Gomez, University College London
                     <c.tur@ucl.ac.uk> <carmen.tur.gomez@gmail.com>
                     
                     Francesco Grussu, University College London
		     <f.grussu@ucl.ac.uk> <francegrussu@gmail.com>
             
                     # Code released under BSD Two-Clause license
                     #
                     # Copyright (c) 2021-2022 University College London. 
                     # All rights reserved.
                     #
                     # Redistribution and use in source and binary forms, with or without modification, are permitted provided that 
                     # the following conditions are met:
                     # 
                     # 1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
                     # 2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer 
                     #    in the documentation and/or other materials provided with the distribution.
                     # 
                     # THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, 
                     # BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. 
                     # IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, 
                     # OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, 
                     # OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, 
                     # OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF 
                     # THE POSSIBILITY OF SUCH DAMAGE.
                     # 
                     # The views and conclusions contained in the software and documentation are those
                     # of the authors and should not be interpreted as representing official policies,
                     # either expressed or implied, of the FreeBSD Project.

	'''
	
	# Get mask dimensions
	mask_size = mask_data.shape
	mask_size = np.array(mask_size)

	# Get voxel positions using the following coordinate system:
	# x = resolution_i*(i - i_ref)
	# y = resolution_j*(j - j_ref)
	# z = resolution_k*(k - k_ref)
	#
	# (i,j,k) are voxel positions in the voxel grid [dimensionless]
	# (resolution_i,resolution_j,resolution_k) are corresponding resolutions [mm]
	# (i_ref,j_ref,k_ref) is the voxel to be used as reference [dimensionless] (example: center-of-mass of motor cortex)
	xpos = np.zeros(mask_size)
	ypos = np.zeros(mask_size)
	zpos = np.zeros(mask_size)

	for ii in range(0,mask_size[0]):
		for jj in range(0,mask_size[1]):
			for kk in range(0,mask_size[2]):
				xpos[ii,jj,kk] = 1.0*mask_resmm[0]*(ii - refvox[0]) 
				ypos[ii,jj,kk] = 1.0*mask_resmm[1]*(jj - refvox[1]) 
				zpos[ii,jj,kk] = 1.0*mask_resmm[2]*(kk - refvox[2])

	# Get positions of all voxels within the lesion mask
	xpos_array = xpos[mask_data==1]   # x-position of voxels belonging to lesions
	ypos_array = ypos[mask_data==1]   # y-position of voxels belonging to lesions
	zpos_array = zpos[mask_data==1]   # z-position of voxels belonging to lesions
	Nvox = xpos_array.size   # Number of voxels
	if(Nvox<=1):
		print('                   WARNING! The lesion mask contains 1 voxel or less and no covariance analysis is possible. Returning NaNs!')
		spacems_indices = np.zeros((1,17))*np.nan
		return spacems_indices 		

	if( (Nvox>1) and (Nvox<10) ):
		print('                   WARNING! The lesion mask contains fewer than 10 voxels. SSPACE-MS analysis will be performed but please, double-check the results!')

	xpos_array = np.reshape(xpos_array,(1,Nvox))       # x-positions: Make sure we store positions and values as row arrays
	ypos_array = np.reshape(ypos_array,(1,Nvox))       # y-positions: Make sure we store positions and values as row arrays
	zpos_array = np.reshape(zpos_array,(1,Nvox))       # z-positions: Make sure we store positions and values as row arrays
	pos_matrix = np.concatenate((xpos_array,ypos_array,zpos_array),axis=0)   # Concatenate x, y and z to have a matrix with rows = variables (x, y, and z) and columns = observations (voxel1, voxel2, etc)

	### Get lesion COM
	pos_com = np.mean(pos_matrix,axis=1)

	#### Study lesion distribution
	# Calculate covariance matrix of lesional voxel positions
	pos_cov = np.cov(pos_matrix)

	# Extract eigenvalues of covariance matrix
	pos_eigval, pos_eigvec = np.linalg.eig(pos_cov)   # Get eigenvalues and eigenvectors
	pos_eigval_sorted = np.array(np.sort(pos_eigval))   # Sort eigenvalues such that pos_eigval_sorted[0] <= pos_eigval_sorted[1] <= pos_eigval_sorted[2]
	
	# Get SSPACE-MS metrics
	axial_cov = pos_eigval_sorted[2]                                                       # Axial covariance (mm^2)
	second_eig = pos_eigval_sorted[1]                                                      # Second eigenvalue (mm^2)
	third_eig = pos_eigval_sorted[0]                                                       # Third eigenvalue (mm^2)
	radial_cov = 0.5*(second_eig + third_eig)                                              # Radial covariance (mm^2)
	mean_cov = (pos_eigval_sorted[0] + pos_eigval_sorted[1] + pos_eigval_sorted[2])/3.0    # Mean covariance (mm^2)
	fa_cov = np.sqrt( (pos_eigval_sorted[0] - mean_cov)**2 + (pos_eigval_sorted[1] - mean_cov)**2 + (pos_eigval_sorted[2] - mean_cov)**2 ) / np.sqrt(pos_eigval_sorted[0]*pos_eigval_sorted[0] + pos_eigval_sorted[1]*pos_eigval_sorted[1] + pos_eigval_sorted[2]*pos_eigval_sorted[2])
	fa_cov = (np.sqrt(3.0/2.0))*(fa_cov)   # Covariance fractional anisotropy (dimensionless)
	tot_les_load = Nvox*mask_resmm[0]*mask_resmm[1]*mask_resmm[2]  # Lesion load (mm^3)
	plan_cov = (2.0/3.0)*(second_eig - third_eig)/(mean_cov)       # Planarity (dimensionless)
	sph_cov = third_eig/mean_cov # Sphericity (dimensionless)
	
	# Return properties of lesion mask (COM and Covariance properties)
	spacems_indices = np.zeros((1,18))
	spacems_indices[0,:] = np.array([axial_cov,second_eig,third_eig,radial_cov,mean_cov,fa_cov,plan_cov,sph_cov,tot_les_load,pos_cov[0,0],pos_cov[1,1],pos_cov[2,2],pos_cov[0,1],pos_cov[0,2],pos_cov[1,2],pos_com[0],pos_com[1],pos_com[2]])
	return spacems_indices




# Function for SSPACE-MS analysis
def runSSPACEMS(input_mask,output_base,input_mask_fu=[],ref_anat=[]):
	''' Perform SSPACE-MS analysis on a binary lesion mask; longitudinal followup masks can also be provided 

	    INTERFACE
	    runSSPACEMS(input_mask,output_base,input_mask_fu=[],ref_anat=[])
	    
	    PARAMETERS
	    - input_mask: input NIFTI file containing the lesion mask (1 in voxels with lesions; 0 otherwise)
	    - output_base: base name of output files, to which the following suffixes will be added: 
                           *_lesion_labels.nii, *_lesion_metrics.csv, *_global_metrics.csv. A description of each of 
                           these is provided below.

                           1) *_lesion_labels.nii is a NIFTI file storing a unique label in each lesion 

                           2) *_lesion_metrics.csv is a CSV file storing SSPACE-MS metrics of each lesion: 
                              - lesid or lesion label; 
                              - l1 aka first eigenvalue aka axcov (axial covariance) (mm^2); 
                              - l2 aka second eigenvalue (mm^2); 
                              - l3 aka third eigenvalue (mm^2); 
                              - radcov aka radial covariance (mm^2); 
                              - meancov aka mean covariance (mm^2); 
                              - covfa aka covariance fractional anisotropy (ranging in [0; 1], with 1 maxiumum anisotropy);
                              - covpl aka covariance planarity (ranging in [0; 1], with 1 maximum planarity); 
                              - covsph aka covariance sphericity (ranging in [0; 1], with 1 maximum sphericity); 
                              - lesvol aka lesion volume (mm^3); 
                              - sxx,syy,szz,sxy,sxz,syz, aka elements of the full covariance matrix (mm^2)
                              (note that x,y,z are image space coordinates); 
                              - comx, comy, comz aka (centre of mass of lesion with respect
			       to anatomical reference if passed, otherwise with respect to the origin of the voxel coordinates; 
                              each of comx, comy, comz is in mm)   

                           3) *_global_metrics.csv is a CSV file storing SSPACE-MS metrics of the global lesion mask: 
                              - l1 aka first eigenvalue aka axcov of axial covariance (mm^2); 
			       - l2 aka second eigenvalue (mm^2); 
			       - l3 aka third eigenvalue (mm^2); 
                              - radcov aka radial covariance (mm^2); 
			       - meancov aka mean covariance (mm^2); 
                              - covfa aka covariance fractional anisotropy (ranging in [0; 1], with 1 maxiumum anisotropy); 
                              - covpl aka covariance planarity (ranging in [0; 1], with 1 maximum planarity);
                              - covsph aka covariance sphericity (ranging in [0; 1], with 1 maximum sphericity);
                              - lesld aka total lesion load (mm^3); 
                              - sxx,syy,szz,sxy,sxz,syz, aka elements of the full covariance matrix (mm^2))
                              (note that x,y,z are image space coordinates); 
                              - comx, comy, comz aka (centre of mass of overall lesion mask with respect
			       to anatomical reference if passed, otherwise with respect to the origin of the voxel coordinates; 
                              each of comx, comy, comz is in mm)   

               

	    - input_mask_fu:  list storing paths to followup NIFTI files storing the followup lesion masks.
                              These followup NIFTI masks should store 1 in voxels with lesions; 0 otherwise.
                              Also, they need to be co-registered to the mask at baseline (input input_mask). 
                              When followups are provided, output files will be produced for each followup point.
                              These will be 
				      *_lesion_labels_tp2.nii, *_lesion_metrics_tp2.csv, *_global_metrics_tp2.csv for timepoint 2;
				      *_lesion_labels_tp3.nii, *_lesion_metrics_tp3.csv, *_global_metrics_tp3.csvfor timepoint 3;
			       Note that followup outputs tp2 will refer to input_mask_followup[0],
                              tp3 will refer to input_mask_followup[1], and so on.   


	    - ref_anat:       path of an optional binary NIFTI file containing the segmentation of a structure to be used as anatomical
			       reference to calculate the center of mass of each lesion. 
			       As an example, this reference file could contain a segmentation of the motor cortex


	    DEPENDENCIES		
		nibabel, numpy, skimage


	    Authors: Carmen Tur Gomez, University College London
                     <c.tur@ucl.ac.uk> <carmen.tur.gomez@gmail.com>
                     
                     Francesco Grussu, University College London
		     <f.grussu@ucl.ac.uk> <francegrussu@gmail.com>
             
                     # Code released under BSD Two-Clause license
                     #
                     # Copyright (c) 2021-2022 University College London. 
                     # All rights reserved.
                     #
                     # Redistribution and use in source and binary forms, with or without modification, are permitted provided that 
                     # the following conditions are met:
                     # 
                     # 1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
                     # 2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer 
                     #    in the documentation and/or other materials provided with the distribution.
                     # 
                     # THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, 
                     # BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. 
                     # IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, 
                     # OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, 
                     # OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, 
                     # OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF 
                     # THE POSSIBILITY OF SUCH DAMAGE.
                     # 
                     # The views and conclusions contained in the software and documentation are those
                     # of the authors and should not be interpreted as representing official policies,
                     # either expressed or implied, of the FreeBSD Project.

	'''

	### Load NIFTI mask
	try:
		mask_obj = nib.load(input_mask)
	except:
		print('')
		print('             ERROR: the 3D lesion mask {} does not exist or is not in NIFTI format. Exiting with 1.'.format(input_mask))
		print('')
		sys.exit(1)		
	mask_data = mask_obj.get_fdata()

	### Make sure the mask is 3D
	mask_header = mask_obj.header
	mask_affine = mask_header.get_best_affine()
	mask_dims = mask_obj.shape
	mask_size = mask_data.shape
	mask_size = np.array(mask_size)
	if mask_size.size!=3:
		print('')
		print('             ERROR: the 3D lesion mask {} is not a 3D NIFTI. Exiting with 1.'.format(input_mask))	 
		print('')
		sys.exit(1)

	### Get image resolution in mm
	mask_resmm = mask_header.get_zooms()
	mask_resmm = np.array(mask_resmm)

	### Check that resolution is roughly isotropic (1% tolerance)
	if( ( (mask_resmm[0]/mask_resmm[1]) > 1.01 ) or ( (mask_resmm[0]/mask_resmm[1]) < 0.99 ) ):
		print('')
		print('             ERROR: the resolution of the 3D lesion mask {} is not isotropic and SSPACE-MS cannot be run.'.format(input_mask))
		print('             Exiting with 1.')
		print('') 
		print('')
		sys.exit(1)
	if( ( (mask_resmm[0]/mask_resmm[2]) > 1.01 ) or ( (mask_resmm[0]/mask_resmm[2]) < 0.99 ) ):
		print('')
		print('             ERROR: the resolution of the 3D lesion mask {} is not isotropic and SSPACE-MS cannot be run.'.format(input_mask))
		print('             Exiting with 1.')
		print('')
		sys.exit(1)

	### Convert mask to binary
	mask_data = np.array(mask_data)
	mask_data[mask_data>0] = 1.0
	mask_data[mask_data<=0] = 0.0
	les_vox = mask_data[mask_data==1]
	les_vox_size = les_vox.size
	if(les_vox_size==0):
		print('             WARNING! The lesion mask is empty (no voxels set to 1) - output metrics will be all NaNs.')

	### Get anatomical reference for Centre-Of-Mass (COM) computation
	if not ref_anat:
		myref_com = [0.0,0.0,0.0] # No reference provided: use origin of voxel grid as reference
	else:
		# Load NIFTI
		try:
			myrefnifti = nib.load(ref_anat)
		except:
			print('')
			print('             ERROR: the 3D anatomical reference {} does not exist or is not in NIFTI format. Exiting with 1.'.format(ref_anat))
			print('')
			sys.exit(1)

		myrefdata = myrefnifti.get_fdata()
		shape_myrefdata = myrefdata.shape		
		header_myrefdata = myrefnifti.header
		affine_myrefdata = header_myrefdata.get_best_affine()			
		if( (np.sum(affine_myrefdata==mask_affine)!=16) or (shape_myrefdata[0]!=mask_dims[0]) or (shape_myrefdata[1]!=mask_dims[1]) or (shape_myrefdata[2]!=mask_dims[2]) ):
			print('             ERROR: the geometry of the anatomical reference {} does not match that of the input mask {}. Exiting with 1.'.format(ref_anat,input_mask))					 
			print('')
			sys.exit(1)

		# Binarise
		myrefdata[myrefdata<=0] = 0.0
		myrefdata[myrefdata>0] = 1.0

		# Extract COM of anatomical reference
		iipos = np.zeros(shape_myrefdata)
		jjpos = np.zeros(shape_myrefdata)
		kkpos = np.zeros(shape_myrefdata)

		for ii in range(0,shape_myrefdata[0]):
			for jj in range(0,shape_myrefdata[1]):
				for kk in range(0,shape_myrefdata[2]):
					iipos[ii,jj,kk] = 1.0*ii 
					jjpos[ii,jj,kk] = 1.0*jj 
					kkpos[ii,jj,kk] = 1.0*kk

		# Get positions of all voxels within the lesion mask
		iipos_array = iipos[myrefdata==1]   # ii-position of voxels belonging to anatomical reference
		jjpos_array = jjpos[myrefdata==1]   # jj-position of voxels belonging to anatomical reference
		kkpos_array = kkpos[myrefdata==1]   # kk-position of voxels belonging to anatomical reference

		# Get reference COM
		myref_com = [np.mean(iipos_array),np.mean(jjpos_array),np.mean(kkpos_array)] # COM of anatomical structure			



	### Load lesion mask at followup time point if provided as optional inputs
	if not (not input_mask_fu):
		Nfu = len(input_mask_fu) # Number of followups
		mask_data_fu_list = []   # List to store the different followup masks
		mask_resmm_fu_list = []  # List to store image resolution

		## Loop through different follow up time points
		for ff in range(0,Nfu):
        
			## Load NIFTI mask at different followups
			try:
				mask_obj_fu = nib.load(input_mask_fu[ff])
			except:
				print('')
				print('             ERROR: the lesion mask at time point {} ({}) does not exist or is not in NIFTI format. Exiting with 1.'.format(ff+2,input_mask_fu[ff]))
				print('')
				sys.exit(1)
			mask_data_fu = mask_obj_fu.get_fdata()
			mask_dims_fu = mask_obj_fu.shape		
			mask_header_fu = mask_obj_fu.header
			mask_affine_fu = mask_header_fu.get_best_affine()			
			if( (np.sum(mask_affine_fu==mask_affine)!=16) or (mask_dims_fu[0]!=mask_dims[0]) or (mask_dims_fu[1]!=mask_dims[1]) or (mask_dims_fu[2]!=mask_dims[2]) ):
				print('             ERROR: the geometry of the mask at time point {} ({}) does not match that of the mask at time point 1 ({}). Exiting with 1.'.format(ff+2,input_mask_fu[ff],input_mask))					 
				print('')
				sys.exit(1)

			### Get image resolution in mm
			mask_resmm_fu = mask_header_fu.get_zooms()
			mask_resmm_fu = np.array(mask_resmm_fu)

			### Convert mask to binary
			mask_data_fu = np.array(mask_data_fu)
			mask_data_fu[mask_data_fu>0] = 1.0
			mask_data_fu[mask_data_fu<=0] = 0.0
			les_vox_fu = mask_data_fu[mask_data_fu==1]
			les_vox_size_fu = les_vox_fu.size
			if(les_vox_size_fu==0):
				print('             WARNING! The lesion mask at time point {} (file {}) is empty (no voxels set to 1) - it looks like all lesions have disappeared!'.format(ff+2,input_mask_fu[ff]))

			### Store 3D matrix with data at followup and release some memory
			mask_data_fu_list.append(mask_data_fu)
			mask_resmm_fu_list.append(mask_resmm_fu)
			del mask_data_fu

 
	#### Get global SSPACE-MS metrics for the input mask
	print('')
	print('SSPACE-MS analysis of overall lesion mask')
	print('')
	print('')
	global_metrics = getSSPACEMSmetrics(mask_data,mask_resmm,myref_com)
	np.savetxt('{}_global_metrics.csv'.format(output_base), global_metrics, fmt='%.6f', delimiter=',', header='axcov,l2,l3,radcov,meancov,covfa,covpl,covsph,lesld,sxx,syy,szz,sxy,sxz,syz,comx,comy,comz',comments='')

	#### Get global SSPACE-MS metrics for the longitudinal input mask (if provided)
	if not (not input_mask_fu):
		for ff in range(0,Nfu):
			print('')
			print('SSPACE-MS analysis of overall lesion mask at followup time point {}'.format(ff+2))
			print('')
			print('')	
			global_metrics_tp2 = getSSPACEMSmetrics(mask_data_fu_list[ff],mask_resmm_fu_list[ff],myref_com)
			np.savetxt('{}_global_metrics_tp{}.csv'.format(output_base,ff+2), global_metrics_tp2, fmt='%.6f', delimiter=',', header='axcov,l2,l3,radcov,meancov,covfa,covpl,covsph,lesld,sxx,syy,szz,sxy,sxz,syz,comx,comy,comz',comments='')



	#### Approximate each lesion to an ellipsoid and study its shape
	print('')
	print('SSPACE-MS analysis of individual lesions                              ')
	print('')
	print('')

	### Extract connected components in lesion mask
	overall_mask = np.copy(mask_data) 
	if not (not input_mask_fu):
		# Get an overall mask by merging multiple followup time points in one joint mask, if these are provided
		for ff in range(0,Nfu):
			overall_mask = overall_mask + mask_data_fu_list[ff]
		overall_mask[overall_mask>1] = 1

	overall_mask_labelled = measure.label(overall_mask)  # Find connected components
	lesion_id = np.unique(overall_mask_labelled)   # Extract lesion labels
	lesion_id = lesion_id[lesion_id>0]             # Remove background
	Nles = lesion_id.size                          # Extract number of lesions
	mask_data_labelled = overall_mask_labelled*mask_data    # Labels for first lesion mask
	lesion_metrics = np.zeros((Nles,19))           # Matrix to store lesion information
	if not (not input_mask_fu):
		mask_data_fu_label_list = []  # List to store lesion labels at followup time points
		lesion_metrics_fu_list = []   # List to store lesion metrics at followup time points
		for ff in range(0,Nfu):
			mask_data_fu_label_list.append(overall_mask_labelled*mask_data_fu_list[ff]) # Labels for followup time points  

		# Release some memory
		del mask_data_fu_list		

	### Study each lesion
	if(Nles==0):
		print('                   WARNING! There are no lesions!')

	# Allocate memory to store lesion-wise SSPACE-MS metrics at followup time points
	if not (not input_mask_fu):
		lesion_metrics_fu = np.zeros((Nfu,Nles,19))

	for qq in range(0,Nles):
		
		print('  ... lesion {} out of {}'.format(qq+1,Nles))
		
		# Get a mask for the qq-th lesion only
		lesqq_mask = np.copy(mask_data_labelled)
		lesqq_mask[lesqq_mask!=lesion_id[qq]] = 0
		lesqq_mask[lesqq_mask>0] = 1

		# Get SSPACE-MS metrics
		print('                   * baseline lesion mask')
		local_metrics_qq = getSSPACEMSmetrics(lesqq_mask,mask_resmm,myref_com)
		lesion_metrics[qq,0] = lesion_id[qq]
		for pp in range(1,19):
			lesion_metrics[qq,pp] = local_metrics_qq[0,pp-1]

		## Study lesions also in the optional longitudinal followup masks (if provided)
		if not (not input_mask_fu):
			for ff in range(0,Nfu):

				print('                   * followup time point {}'.format(ff+2))

				# Get a mask with current lesion only
				lesqq_mask_tp2 = np.copy(mask_data_fu_label_list[ff])
				lesqq_mask_tp2[lesqq_mask_tp2!=lesion_id[qq]] = 0
				lesqq_mask_tp2[lesqq_mask_tp2>0] = 1

				# Get SSPACE-MS metrics
				local_metrics_qq_tp2 = getSSPACEMSmetrics(lesqq_mask_tp2,mask_resmm_fu_list[ff],myref_com)
				lesion_metrics_fu[ff,qq,0] = lesion_id[qq]
				for pp in range(1,19):
					lesion_metrics_fu[ff,qq,pp] = local_metrics_qq_tp2[0,pp-1]
				
				# Store lesion-wise metrics for current followup in the overall lesion list for followup masks
				lesion_metrics_fu_list.append(np.ndarray.squeeze(lesion_metrics_fu[ff,:,:]))  

		
	print('')	
			
	### Save lesion-wise metrics
	# Save lesion-wise SSPACE-MS metrics to an output text file
	np.savetxt('{}_lesion_metrics.csv'.format(output_base), lesion_metrics, fmt='%.6f', delimiter=',', header='lesid,axcov,l2,l3,radcov,meancov,covfa,covpl,covsph,lesvol,sxx,syy,szz,sxy,sxz,syz,comx,comy,comz',comments='')
	
	# Save lesion labels as a NIFTI file
	buffer_header = mask_obj.header
	buffer_header.set_data_dtype('float32')   # Make sure we save output data as float32, even if input header indicates a different data type
	label_obj = nib.Nifti1Image(mask_data_labelled,mask_obj.affine,buffer_header)
	nib.save(label_obj, '{}_lesion_labels.nii'.format(output_base))

	### Save lesion-wise metrics for the longitudinal followup masks
	if not (not input_mask_fu):

		for ff in range(0,Nfu): 
			# Save lesion-wise SSPACE-MS metrics to an output text file
			np.savetxt('{}_lesion_metrics_tp{}.csv'.format(output_base,ff+2), lesion_metrics_fu_list[ff], fmt='%.6f', delimiter=',', header='lesid,axcov,l2,l3,radcov,meancov,covfa,covpl,covsph,lesvol,sxx,syy,szz,sxy,sxz,syz,comx,comy,comz',comments='')
			
			# Save lesion labels as a NIFTI file
			buffer_header = mask_obj.header
			buffer_header.set_data_dtype('float32')   # Make sure we save output data as float32, even if input header indicates a different data type
			label_obj_tp2 = nib.Nifti1Image(mask_data_fu_label_list[ff],mask_obj.affine,buffer_header)
			nib.save(label_obj_tp2, '{}_lesion_labels_tp{}.nii'.format(output_base,ff+2))




# Run the module as a script when required
if __name__ == "__main__":

	### Print help and parse arguments
	parser = argparse.ArgumentParser(description='Perform spatial analysis of lesions following the SSPACE-MS framework. Dependencies: nibabel, numpy, skimage. Authors: Carmen Tur <c.tur@ucl.ac.uk> <carmen.tur.gomez@gmail.com> and Francesco Grussu <f.grussu@ucl.ac.uk> <francegrussu@gmail.com>. # Code released under BSD Two-Clause license # Copyright (c) 2021-2022 University College London. # All rights reserved. ## Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met: ## 1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer. # 2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution. ## THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE. ## The views and conclusions contained in the software and documentation are those of the authors and should not be interpreted as representing official policies, either expressed or implied, of the FreeBSD Project.')
	parser.add_argument('mask_nifti', help='input NIFTI file containing the lesion mask (1 in voxels flagging lesions; 0 otherwise)')
	parser.add_argument('out_basename', help='base name of output files, to which the following suffixes will be added: *_lesion_labels.nii, *_lesion_metrics.csv, *_global_metrics.csv. *** FILE *_lesion_labels.nii is a NIFTI file storing a unique label in each lesion; *** FILE *_lesion_metrics.csv is a CSV file storing SSPACE-MS metrics of each lesion: lesid or lesion label; axcov of axial covariance (mm^2); l2 or second eigenvalue (mm^2); l3 or third eigenvalue (mm^2); radcov or radial covariance (mm^2); meancov or mean covariance (mm^2); covfa or covariance fractional anisotropy (ranging in [0; 1], with 1 maxiumum anisotropy); covpl or covariance planarity (ranging in [0; 1], with 1 maximum planarity); lesvol or lesion volume (mm^3); sxx,syy,szz,sxy,sxz,syz, elements of the full covariance matrix (mm^2); comx,comy,comz, centre-of-mass of individual lesion with respect to anatomical reference (mm). *** FILE *_global_metrics.csv is a CSV file storing SSPACE-MS metrics of the global lesion mask: axcov of axial covariance (mm^2); l2 or second eigenvalue (mm^2); l3 or third eigenvalue (mm^2); radcov or radial covariance (mm^2); meancov or mean covariance (mm^2); covfa or covariance fractional anisotropy (ranging in [0; 1], with 1 maxiumum anisotropy); covpl or covariance planarity (ranging in [0; 1], with 1 maximum planarity); covsph or covariance sphericity (ranging in [0; 1], with 1 maximum sphericity); lesld or total lesion load (mm^3); sxx,syy,szz,sxy,sxz,syz, elements of the full covariance matrix (mm^2)); comx,comy,comz, centre-of-mass of lesion mask with respect to anatomical reference (mm).')
	parser.add_argument('--fu', metavar='<NIFTI2,NIFTI3,...>', help='list of input NIFTI file containing lesion masks at followups (1 in voxels with lesions; 0 otherwise - note that the followup masks need to be co-registered to that at baseline). If more than one followup is available, these are listed one after the other separating with commas (,). As a consequence, no commas can be contained in the file path. Output files are produced for each followup as well (*_lesion_labels_tp2.nii, *_lesion_metrics_tp2.csv, *_global_metrics_tp2.csv; *_lesion_labels_tp3.nii, *_lesion_metrics_tp3.csv, *_global_metrics_tp3.csv; etc, with tp2 referring to the first file; tp3 to the second file; etc)')
	parser.add_argument('--ref', metavar='<NIFTI>', help='path of an optional binary NIFTI file containing the segmentation of a structure to be used as anatomical reference to calculate the center of mass of each lesion. As an example, this reference file could contain a segmentation of the motor cortex: the centre of mass of each lesion would then be expressed using the centre of mass of the motor cortex as origin of the coordinate system. When this optional input is not provided, the centre of mass of the lesion is expressed using the origin of the voxel grid as reference')
	args = parser.parse_args()

	### Get input arguments
	maskfile = args.mask_nifti
	outname = args.out_basename
	mask2file = args.fu
	refanatfile = args.ref
    
	### Deal with optional arguments
	if isinstance(mask2file, str)==1:
		# Followup masks have been provided
		maskrequest2 = True
		masks_list = (args.fu).split(',')
	else:
		# Followup masks have not been provided
		maskrequest2 = False

	if isinstance(refanatfile, str)==1:
		# Anatomical reference has been provided
		refrequest = True
	else:
		refrequest = False


	### Print a message
	print('')
	print('*****************************************************************************')
	print('                             SSPACE-MS analysis                              ')
	print('*****************************************************************************')
	print('')
	print('Input NIFTI mask: {}'.format([maskfile]))
	print('Output base name: {}'.format([outname]))
	if (maskrequest2==True):
		print('Followup masks: {}'.format(masks_list))
	if (refrequest==True):
		print('Anatomical reference for CoM calculation: {}'.format([refanatfile]))


	### Run SSPACE-MS analysis
	if((maskrequest2==False) and (refrequest==False)):
		print('')
		runSSPACEMS(maskfile,outname)

	elif((maskrequest2==True) and (refrequest==False)):
		print('')
		runSSPACEMS(maskfile,outname,input_mask_fu=masks_list)

	elif((maskrequest2==False) and (refrequest==True)):
		print('')
		runSSPACEMS(maskfile,outname,ref_anat=refanatfile)

	elif((maskrequest2==True) and (refrequest==True)):
		print('')
		runSSPACEMS(maskfile,outname,input_mask_fu=masks_list,ref_anat=refanatfile)



	### Done
	print('')
	print('')
	print('Processing completed.')
	print('')
	sys.exit(0)


