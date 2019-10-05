# SPACE-MS analysis
#
# Authors:          Carmen Tur Gomez, University College London
#                  <c.tur@ucl.ac.uk> <carmen.tur.gomez@gmail.com>
#                   Francesco Grussu, University College London
#		   <f.grussu@ucl.ac.uk> <francegrussu@gmail.com>
#
# Code released under BSD Two-Clause license
#
# Copyright (c) 2018 University College London. 
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

# Function for SPACE-MS analysis
def runSPACEMS(input_mask,output_base):
	''' Perform SPACE-MS analysis on a binary lesion mask  

	    INTERFACE
	    runSPACEMS(input_mask,output_base)
	    
	    PARAMETERS
	    - input_mask: input NIFTI file containing the lesion mask (1 in voxels flagging lesions; 0 otherwise)
	    - output_base: base name of output files, to which the following suffixes will be added: 
                           *_lesion_labels.nii, *_lesion_metrics.csv, *_global_metrics.csv. 

                           1) *_lesion_labels.nii is a NIFTI file storing a unique label in each lesion 

                           2) *_lesion_metrics.csv is a CSV file storing SPACE-MS metrics of each lesion: 
                              lesid or lesion label; axcov of axial covariance (mm^2); l2 or second eigenvalue (mm^2); 
                              l3 or third eigenvalue (mm^2); radcov or radial covariance (mm^2); meancov or mean covariance (mm^2); 
                              covfa or covariance fractional anisotropy (ranging in [0; 1], with 1 maxiumum anisotropy); 
                              covpl or covariance planarity (ranging in [-1; 1], with -1 maximum planarity); 
                              lesvol or lesion volume (mm^3); sxx,syy,szz,sxy,sxz,syz, elements of the full covariance matrix (mm^2) 

                           3) *_global_metrics.csv is a CSV file storing SPACE-MS metrics of the global lesion mask: 
                              axcov of axial covariance (mm^2); l2 or second eigenvalue (mm^2); l3 or third eigenvalue (mm^2); 
                              radcov or radial covariance (mm^2); meancov or mean covariance (mm^2); 
                              covfa or covariance fractional anisotropy (ranging in [0; 1], with 1 maxiumum anisotropy); 
                              covpl or covariance planarity (ranging in [-1; 1], with -1 maximum planarity); 
                              lesld or total lesion load (mm^3); 
                              sxx,syy,szz,sxy,sxz,syz, elements of the full covariance matrix (mm^2)).

	    Authors: Carmen Tur Gomez, University College London
                     <c.tur@ucl.ac.uk> <carmen.tur.gomez@gmail.com>
                     
                     Francesco Grussu, University College London
		     <f.grussu@ucl.ac.uk> <francegrussu@gmail.com>
             
                     # Code released under BSD Two-Clause license
                     #
                     # Copyright (c) 2018 University College London. 
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
		print('ERROR: the 3D lesion mask {} does not exist or is not in NIFTI format. Exiting with 1.'.format(input_mask))
		print('')
		sys.exit(1)		
	mask_data = mask_obj.get_data()

	### Make sure the mask is 3D
	mask_header = mask_obj.header
	mask_affine = mask_header.get_best_affine()
	mask_dims = mask_obj.get_shape()
	mask_size = mask_data.shape
	mask_size = np.array(mask_size)
	if mask_size.size!=3:
		print('')
		print('ERROR: the 3D lesion mask {} is not a 3D NIFTI. Exiting with 1.'.format(input_mask))	 
		print('')
		sys.exit(1)

	### Get image resolution in mm
	mask_resmm = mask_header.get_zooms()
	mask_resmm = np.array(mask_resmm)

	### Check that resolution is isotropic
	if mask_resmm[0]!=mask_resmm[1]:
		print('')
		print('ERROR: the resolution of the 3D lesion mask {} is not isotropic - metrics would therefore be biased.'.format(input_mask))
		print('Exiting with 1.')
		print('') 
		print('')
		sys.exit(1)
	if mask_resmm[0]!=mask_resmm[2]:
		print('')
		print('ERROR: the resolution of the 3D lesion mask {} is not isotropic - metrics would therefore be biased.'.format(input_mask))
		print('Exiting with 1.')
		print('')
		sys.exit(1)

	### Convert mask to binary
	mask_data = np.array(mask_data)
	mask_data[mask_data>0] = 1.0
	mask_data[mask_data<=0] = 0.0
	les_vox = mask_data[mask_data==1]
	les_vox_size = les_vox.size
	if(les_vox_size==0):
		print('ERROR! The lesion mask is empty (no voxels set to 1). Exiting with 1, no output files saved!')
		print('')
		sys.exit(1)

	# Get voxel positions
	xpos = np.zeros(mask_size)
	ypos = np.zeros(mask_size)
	zpos = np.zeros(mask_size)

	for ii in range(0,mask_size[0]):
		for jj in range(0,mask_size[1]):
			for kk in range(0,mask_size[2]):
				xpos[ii,jj,kk] = mask_resmm[0]*ii*1.0
				ypos[ii,jj,kk] = mask_resmm[1]*jj*1.0
				zpos[ii,jj,kk] = mask_resmm[2]*kk*1.0

	# Get positions of all voxels within the lesion mask
	xpos_array = xpos[mask_data==1]   # x-position of voxels belonging to lesions
	ypos_array = ypos[mask_data==1]   # y-position of voxels belonging to lesions
	zpos_array = zpos[mask_data==1]   # z-position of voxels belonging to lesions
	Nvox = xpos_array.size   # Number of voxels
	if(Nvox==1):
		print('ERROR! The lesion mask contains only 1 voxel and no covariance analysis is possible. Exiting with 1, no output files saved!')
		print('')
		sys.exit(1)
	if( (Nvox>1) and (Nvox<10) ):
		print('WARNING! The lesion mask contains less than 10 voxels. SPACE-MS analysis will be performed but please, double-check the results!')
		print('')

	xpos_array = np.reshape(xpos_array,(1,Nvox))       # x-positions: Make sure we store positions and values as row arrays
	ypos_array = np.reshape(ypos_array,(1,Nvox))       # y-positions: Make sure we store positions and values as row arrays
	zpos_array = np.reshape(zpos_array,(1,Nvox))       # z-positions: Make sure we store positions and values as row arrays
	pos_matrix = np.concatenate((xpos_array,ypos_array,zpos_array),axis=0)   # Concatenate x, y and z to have a matrix with rows = variables (x, y, and z) and columns = observations (voxel1, voxel2, etc)

	#### Study the overall distribution of MS damage in the brain

	# Calculate covariance matrix of lesional voxel positions
	pos_cov = np.cov(pos_matrix)

	# Extract eigenvalues of covariance matrix
	pos_eigval, pos_eigvec = np.linalg.eig(pos_cov)   # Get eigenvalues and eigenvectors
	pos_eigval_sorted = np.array(np.sort(pos_eigval))   # Sort eigenvalues such that pos_eigval_sorted[0] <= pos_eigval_sorted[1] <= pos_eigval_sorted[2]
	
	# Get SPACE-MS metrics
	axial_cov = pos_eigval_sorted[2]                                                       # Axial covariance (mm^2)
	second_eig = pos_eigval_sorted[1]                                                      # Second eigenvalue (mm^2)
	third_eig = pos_eigval_sorted[0]                                                       # Third eigenvalue (mm^2)
	radial_cov = 0.5*(second_eig + third_eig)                                              # Radial covariance (mm^2)
	mean_cov = (pos_eigval_sorted[0] + pos_eigval_sorted[1] + pos_eigval_sorted[2])/3.0    # Mean covariance (mm^2)
	fa_cov = np.sqrt( (pos_eigval_sorted[0] - mean_cov)**2 + (pos_eigval_sorted[1] - mean_cov)**2 + (pos_eigval_sorted[2] - mean_cov)**2 ) / np.sqrt(pos_eigval_sorted[0]*pos_eigval_sorted[0] + pos_eigval_sorted[1]*pos_eigval_sorted[1] + pos_eigval_sorted[2]*pos_eigval_sorted[2])
	fa_cov = (np.sqrt(3.0/2.0))*(fa_cov)   # Covariance fractional anisotropy (dimensionless)
	tot_les_load = Nvox*mask_resmm[0]*mask_resmm[1]*mask_resmm[2]  # Lesion load (mm^3)
	pos_cov_dev = np.copy(pos_cov)  # Allocate deviatoric covariance matrix
	pos_cov_dev[0,0] = pos_cov_dev[0,0] - mean_cov   # Calculate deviatoric component x,x
	pos_cov_dev[1,1] = pos_cov_dev[1,1] - mean_cov   # Calculate deviatoric component y,y
	pos_cov_dev[2,2] = pos_cov_dev[2,2] - mean_cov   # Calculate deviatoric component z,z
	pos_cov_dev_norm = np.sqrt( np.trace( np.matmul(pos_cov_dev,np.matrix.transpose(pos_cov_dev)) ) ) # Norm of deviatoric tensor
	if(pos_cov_dev_norm==0):
		plan_cov = np.nan  # No planarity can be calculated for isotropic tensors
		print('WARNING! The planarity of the covariance matrix relateive to the whole lesion mask is note defined and will be saved as nan!')
		print('')
	else:
		plan_cov = (3*np.sqrt(6))*np.linalg.det(pos_cov_dev/pos_cov_dev_norm) # Planarity (dimensionless)
	

	# Print global properties of lesion mask
	print('')
	print('GLOBAL PROPERTIES OF TISSUE DAMAGE                              ')
	print('')
	print('')
	print('     AC or axial covariance (first eigenvalue) = {} mm^2'.format(axial_cov))
	print('     Second eigenvalue = {} mm^2'.format(second_eig))
	print('     Third eigenvalue = {} mm^2'.format(third_eig))
	print('     RC or radial covariance = {} mm^2'.format(radial_cov))
	print('     MC or mean covariance = {} mm^2'.format(mean_cov))
	print('     CFA or covariance fractional anisotropy = {}'.format(fa_cov))
	print('     CP or covariance planarity = {}'.format(plan_cov))
	print('     Total lesion load = {} mm^3'.format(tot_les_load))
	print('')
	print('')
	print('')

	# Save global properties of lesions to text file
	global_metrics = np.zeros((1,14))
	global_metrics[0,:] = np.array([axial_cov,second_eig,third_eig,radial_cov,mean_cov,fa_cov,plan_cov,tot_les_load,pos_cov[0,0],pos_cov[1,1],pos_cov[2,2],pos_cov[0,1],pos_cov[0,2],pos_cov[1,2]])

	np.savetxt('{}_global_metrics.csv'.format(output_base), global_metrics, fmt='%.6f', delimiter=',', header='axcov,l2,l3,radcov,meancov,covfa,covpl,lesld,sxx,syy,szz,sxy,sxz,syz',comments='')



	#### Approximate each lesion to an ellipsoid and study its shape
	print('')
	print('PROPERTIES OF INDIVIDUAL LESIONS                              ')
	print('')
	print('')

	# Extract connected components in lesion mask
	mask_data_labelled = measure.label(mask_data)  # Find connected components
	lesion_id = np.unique(mask_data_labelled)      # Extract lesion labels
	lesion_id = lesion_id[lesion_id>0]             # Remove background
	Nles = lesion_id.size                          # Extract number of lesions
	
	# Create matrices to store lesion information
	lesion_metrics = np.zeros((Nles,15))

	# Study each lesion
	for qq in range(0,Nles):
		
		# Get positions of voxels within current lesion
		xpos_array_qq = xpos[mask_data_labelled==lesion_id[qq]]            # x-position of voxels belonging to current lesion
		ypos_array_qq = ypos[mask_data_labelled==lesion_id[qq]]            # y-position of voxels belonging to current lesion
		zpos_array_qq = zpos[mask_data_labelled==lesion_id[qq]]            # z-position of voxels belonging to current lesion
		Nvox_qq = xpos_array_qq.size    # Number of voxels in current lesion

		# The lesion has only one voxel: provide nan for SPACE-MS metrics		
		if(Nvox_qq==1):
			print('WARNING! Lesion number {} contains only 1 voxel and no covariance analysis is possible. Its SPACE-MS metrics will be set to nan'.format(lesion_id[qq]))
			print('')

			axial_cov_qq = np.nan    # Axial covariance (mm^2)
			second_eig_qq = np.nan   # Second eigenvalue (mm^2)
			third_eig_qq = np.nan    # Third eigenvalue (mm^2)
			radial_cov_qq = np.nan   # Radial covariance (mm^2)
			mean_cov_qq = np.nan     # Mean covariance (mm^2)
			fa_cov_qq = np.nan       # Covariance fractional anisotropy (dimensionless)
			plan_cov_qq = np.nan     # Planarity (dimensionless) 
			les_vol_qq = 1.0*mask_resmm[0]*mask_resmm[1]*mask_resmm[2]  # Volume of current lesion (mm^3)
			pos_cov_qq =  np.ones((3,3))  
			pos_cov_qq[:] = np.nan     # Full covariance matrix

			# Store SPACE-MS metrics of current lesion to a matrix that will be uses to save the output file
			lesion_metrics[qq,:] = np.array([lesion_id[qq],axial_cov_qq,second_eig_qq,third_eig_qq,radial_cov_qq,mean_cov_qq,fa_cov_qq,plan_cov_qq,les_vol_qq,pos_cov_qq[0,0],pos_cov_qq[1,1],pos_cov_qq[2,2],pos_cov_qq[0,1],pos_cov_qq[0,2],pos_cov_qq[1,2]])

		# The lesion has got more than one voxel: calculate SPACE-MS metrics, and print a warning when lesion volume is less than 10 voxels
		else:
			if( (Nvox_qq>1) and (Nvox_qq<10) ):
				print('WARNING! Lesion number {} contains less than 10 voxels. SPACE-MS analysis will be performed but please, double-check the results!'.format(lesion_id[qq]))
				print('')

			# Prepare positions for covariance analysis
			xpos_array_qq = np.reshape(xpos_array_qq,(1,Nvox_qq))       # x-positions: Make sure we store positions and values as row arrays
			ypos_array_qq = np.reshape(ypos_array_qq,(1,Nvox_qq))       # y-positions: Make sure we store positions and values as row arrays
			zpos_array_qq = np.reshape(zpos_array_qq,(1,Nvox_qq))       # z-positions: Make sure we store positions and values as row arrays
			pos_matrix_qq = np.concatenate((xpos_array_qq,ypos_array_qq,zpos_array_qq),axis=0)   # Concatenate x, y and z to have a matrix with rows = variables (x, y, and z) and columns = observations (voxel1, voxel2, etc)

			# Calculate covariance matrix of lesional voxel positions
			pos_cov_qq = np.cov(pos_matrix_qq)

			# Extract eigenvalues of covariance matrix
			pos_eigval_qq, pos_eigvec_qq = np.linalg.eig(pos_cov_qq)            # Get eigenvalues and eigenvectors
			pos_eigval_sorted_qq = np.array(np.sort(pos_eigval_qq))   # Sort eigenvalues such that pos_eigval_sorted[0] <= pos_eigval_sorted[1] <= pos_eigval_sorted[2]

			# Calculate SPACE-MS metrics of current lesion
			axial_cov_qq = pos_eigval_sorted_qq[2]                                                             # Axial covariance (mm^2)
			second_eig_qq = pos_eigval_sorted_qq[1]                                                            # Second eigenvalue (mm^2)
			third_eig_qq = pos_eigval_sorted_qq[0]                                                             # Third eigenvalue (mm^2)
			radial_cov_qq = 0.5*(second_eig_qq + third_eig_qq)                                                 # Radial covariance (mm^2)
			mean_cov_qq = (pos_eigval_sorted_qq[0] + pos_eigval_sorted_qq[1] + pos_eigval_sorted_qq[2])/3.0    # Mean covariance (mm^2)
			fa_cov_qq = np.sqrt( (pos_eigval_sorted_qq[0] - mean_cov_qq)**2 + (pos_eigval_sorted_qq[1] - mean_cov_qq)**2 + (pos_eigval_sorted_qq[2] - mean_cov_qq)**2 ) / np.sqrt(pos_eigval_sorted_qq[0]*pos_eigval_sorted_qq[0] + pos_eigval_sorted_qq[1]*pos_eigval_sorted_qq[1] + pos_eigval_sorted_qq[2]*pos_eigval_sorted_qq[2])
			fa_cov_qq = (np.sqrt(3.0/2.0))*(fa_cov_qq)   # Covariance fractional anisotropy (dimensionless)
			les_vol_qq = Nvox_qq*mask_resmm[0]*mask_resmm[1]*mask_resmm[2]  # Volume of current lesion (mm^3)

			pos_cov_dev_qq = np.copy(pos_cov_qq)  # Allocate deviatoric covariance matrix
			pos_cov_dev_qq[0,0] = pos_cov_dev_qq[0,0] - mean_cov_qq   # Calculate deviatoric component x,x
			pos_cov_dev_qq[1,1] = pos_cov_dev_qq[1,1] - mean_cov_qq   # Calculate deviatoric component y,y
			pos_cov_dev_qq[2,2] = pos_cov_dev_qq[2,2] - mean_cov_qq   # Calculate deviatoric component z,z
			pos_cov_dev_norm_qq = np.sqrt( np.trace( np.matmul(pos_cov_dev_qq,np.matrix.transpose(pos_cov_dev_qq)) ) ) # Norm of deviatoric tensor
			if(pos_cov_dev_norm_qq==0):
				plan_cov_qq = np.nan  # No planarity can be calculated for isotropic tensors
				print('WARNING! The planarity of lesion {} is not defined and will be saved as nan'.format(lesion_id[qq]))
				print('')
			else:
				plan_cov_qq = (3*np.sqrt(6))*np.linalg.det(pos_cov_dev_qq/pos_cov_dev_norm_qq)   # Planarity (dimensionless) 
		
			# Store SPACE-MS metrics of current lesion to a matrix that will be uses to save the output file
			lesion_metrics[qq,:] = np.array([lesion_id[qq],axial_cov_qq,second_eig_qq,third_eig_qq,radial_cov_qq,mean_cov_qq,fa_cov_qq,plan_cov_qq,les_vol_qq,pos_cov_qq[0,0],pos_cov_qq[1,1],pos_cov_qq[2,2],pos_cov_qq[0,1],pos_cov_qq[0,2],pos_cov_qq[1,2]])

			# Print properties of current lesion
			print('     Lesion ID = {}'.format(lesion_id[qq]))
			print('')
			print('             AC or axial covariance (first eigenvalue) = {} mm^2'.format(axial_cov_qq))
			print('             Second eigenvalue = {} mm^2'.format(second_eig_qq))
			print('             Third eigenvalue = {} mm^2'.format(third_eig_qq))
			print('             RC or radial covariance = {} mm^2'.format(radial_cov_qq))
			print('             MC or mean covariance = {} mm^2'.format(mean_cov_qq))
			print('             CFA or covariance fractional anisotropy = {}'.format(fa_cov_qq))
			print('             CP or covariance planarity = {}'.format(plan_cov_qq))
			print('             Lesion volume = {} mm^3'.format(les_vol_qq))
			print('')
			print('')
			print('')
			

	# Save lesion-wise SPACE-MS metrics to an output text file
	np.savetxt('{}_lesion_metrics.csv'.format(output_base), lesion_metrics, fmt='%.6f', delimiter=',', header='lesid,axcov,l2,l3,radcov,meancov,covfa,covpl,lesvol,sxx,syy,szz,sxy,sxz,syz',comments='')
	
	# Save lesion labels as a NIFTI file
	buffer_header = mask_obj.header
	buffer_header.set_data_dtype('float64')   # Make sure we save output data as float64, even if input header indicates a different data type
	label_obj = nib.Nifti1Image(mask_data_labelled,mask_obj.affine,buffer_header)
	nib.save(label_obj, '{}_lesion_labels.nii'.format(output_base))






# Run the module as a script when required
if __name__ == "__main__":

	### Print help and parse arguments
	parser = argparse.ArgumentParser(description='Perform spatial analysis of lesions following the SPACE-MS framework. Authors: Carmen Tur <c.tur@ucl.ac.uk> <carmen.tur.gomez@gmail.com> and Francesco Grussu <f.grussu@ucl.ac.uk> <francegrussu@gmail.com>. # Code released under BSD Two-Clause license # Copyright (c) 2018 University College London. # All rights reserved. ## Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met: ## 1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer. # 2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution. ## THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE. ## The views and conclusions contained in the software and documentation are those of the authors and should not be interpreted as representing official policies, either expressed or implied, of the FreeBSD Project.')
	parser.add_argument('mask_nifti', help='input NIFTI file containing the lesion mask (1 in voxels flagging lesions; 0 otherwise)')
	parser.add_argument('out_basename', help='base name of output files, to which the following suffixes will be added: *_lesion_labels.nii, *_lesion_metrics.csv, *_global_metrics.csv. *_lesion_labels.nii is a NIFTI file storing a unique label in each lesion; *_lesion_metrics.csv is a CSV file storing SPACE-MS metrics of each lesion: lesid or lesion label; axcov of axial covariance (mm^2); l2 or second eigenvalue (mm^2); l3 or third eigenvalue (mm^2); radcov or radial covariance (mm^2); meancov or mean covariance (mm^2); covfa or covariance fractional anisotropy (ranging in [0; 1], with 1 maxiumum anisotropy); covpl or covariance planarity (ranging in [-1; 1], with -1 maximum planarity); lesvol or lesion volume (mm^3); sxx,syy,szz,sxy,sxz,syz, elements of the full covariance matrix (mm^2); *_global_metrics.csv is a CSV file storing SPACE-MS metrics of the global lesion mask: axcov of axial covariance (mm^2); l2 or second eigenvalue (mm^2); l3 or third eigenvalue (mm^2); radcov or radial covariance (mm^2); meancov or mean covariance (mm^2); covfa or covariance fractional anisotropy (ranging in [0; 1], with 1 maxiumum anisotropy); covpl or covariance planarity (ranging in [-1; 1], with -1 maximum planarity); lesld or total lesion load (mm^3); sxx,syy,szz,sxy,sxz,syz, elements of the full covariance matrix (mm^2)).')
	args = parser.parse_args()

	### Get input arguments
	maskfile = args.mask_nifti
	outname = args.out_basename

	### Print a message
	print('')
	print('*****************************************************************************')
	print('                            SPACE-MS analysis                          ')
	print('*****************************************************************************')
	print('')
	print('Input NIFTI file: {}'.format(maskfile))
	print('Output files: {}_global_metrics.csv {}_lesion_metrics.csv {}_lesion_labels.nii'.format(outname,outname,outname))
	print('')

	### Run SPACE-MS analysis
	runSPACEMS(maskfile,outname)
	
	### Done
	print('Processing completed.')
	print('')
	sys.exit(0)


