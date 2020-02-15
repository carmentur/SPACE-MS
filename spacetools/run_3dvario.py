# VARIO-MS analysis
#
# Authors:          Carmen Tur Gomez, University College London
#                  <c.tur@ucl.ac.uk> <carmen.tur.gomez@gmail.com>
#                   Francesco Grussu, University College London
#		   <f.grussu@ucl.ac.uk> <francegrussu@gmail.com>
#
# Code released under BSD Two-Clause license
#
# Copyright (c) 2020 University College London. 
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
from pykrige.ok3d import OrdinaryKriging3D 


# Function for VARIO-MS analysis based on variogram fitting
def runVARIOMS(input_nifti,mask_nifti,output_base,lag_number,boot_number,vox_number):
	''' Perform exponential variogram analysis on a 3D metric within a 3D mask  

	    INTERFACE
	    runVARIOMS(input_nifti,mask_nifti,output_base)
	    
	    PARAMETERS
	    - input_nifti: input NIFTI file containing a metric whose spatial variability is of interest
                           (for instance, a lesion mask, storing 1 in voxels flagging lesions; 0 otherwise)

	    - mask_nifti: input NIFTI file storing a binary wask higlighting the anatomical region of interest
                           (for instance, the brain)

	    - output_base: base name of output files, to which the following suffixes will be added: *_semivar.csv, *_lags.csv, *_fit.csv. 
                           
		           *_semivar.csv stores the semivariance across bootstrap iterations (a matrix whose rows are bootstrap 
                           iterations and whose columns show values of semivariance for increasing lag). The units of the elements of the matrix 
                           are those of the input metric (input_nifti) to the power of two. 
                          
                           *_lags.csv stores the lags across bootstrap iterations corresponding to the values of semivariance 
                           reported in *_semivar.csv  (it is again a matrix whose rows are bootstrap iterations and whose columns show 
                           the lags at which semivariances are calculated). The units of the elements of the matrix are [mm]. 
                          
                           *_fit.csv is a CSV file storing the parameters of an exponential fit of the semivariance and lag data. 
                           Rows correspond to different bootstrap iterations, while columns correspond to bootstrap index; 
                           sill (same units as semivariance), range (in [mm]), nugget (same units as semivariance).

	     - lag_number: number of lags for variogram fitting

	     - boot_number: number of bootstrap variogram calculations

	     - vox_number: number of voxels subsampled at random to be used for variogram fitting in each bootstrap calculation

	    Authors: Carmen Tur Gomez, University College London
                     <c.tur@ucl.ac.uk> <carmen.tur.gomez@gmail.com>
                     
                     Francesco Grussu, University College London
		     <f.grussu@ucl.ac.uk> <francegrussu@gmail.com>
             
                     # Code released under BSD Two-Clause license
                     #
                     # Copyright (c) 2020 University College London. 
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

	print('')	
	print('Loading input files')
	print('')

	### Load input NIFTI
	try:
		input_obj = nib.load(input_nifti)
	except:
		print('')
		print('ERROR: the 3D input NIFTI {} does not exist or is not in NIFTI format. Exiting with 1.'.format(input_nifti))
		print('')
		sys.exit(1)		
	input_data = input_obj.get_data()

	### Make sure the input is 3D
	input_header = input_obj.header
	input_affine = input_header.get_best_affine()
	input_dims = input_obj.shape
	input_size = input_data.shape
	input_size = np.array(input_size)
	if input_size.size!=3:
		print('')
		print('ERROR: the 3D input NIFTI {} is not 3D. Exiting with 1.'.format(input_mask))	 
		print('')
		sys.exit(1)

	### Get the region-of-interest mask
	try:
		roi_obj = nib.load(mask_nifti)
	except:
		print('')
		print('ERROR: the 3D mask NIFTI {} does not exist or is not in NIFTI format. Exiting with 1.'.format(mask_nifti))
		print('')
		sys.exit(1)		
	roi_data = roi_obj.get_data()

	### Make sure the ROI NIFTI is 3D
	roi_header = roi_obj.header
	roi_affine = roi_header.get_best_affine()
	roi_dims = roi_obj.shape
	roi_size = roi_data.shape
	roi_size = np.array(roi_size)
	if roi_size.size!=3:
		print('')
		print('ERROR: the 3D mask NIFTI {} is not 3D. Exiting with 1.'.format(mask_nifti))	 
		print('')
		sys.exit(1)

	### Make sure the ROI NIFTI and the input NIFTI have consistent headers
	if ( (np.sum(input_affine==roi_affine)!=16) or (input_dims[0]!=roi_dims[0]) or (input_dims[1]!=roi_dims[1]) or (input_dims[2]!=roi_dims[2]) ):
			print('')
			print('ERROR: the geometry of the mask file {} does not match that of the input data. Exiting with 1.'.format(mask_nifti))					 
			print('')
			sys.exit(1)


	### Make sure numeric inputs are numeric
	try:
		lag_number = float(lag_number)
	except:
		print('')
		print('ERROR: the lag number {} does not appear to be a numeric scalar. Exiting with 1.'.format(lag_number))
		print('')
		sys.exit(1)

	try:
		boot_number = float(boot_number)
	except:
		print('')
		print('ERROR: the bootstrap number {} does not appear to be a numeric scalar. Exiting with 1.'.format(boot_number))
		print('')
		sys.exit(1)

	try:
		vox_number = float(vox_number)
	except:
		print('')
		print('ERROR: the voxel number {} does not appear to be a numeric scalar. Exiting with 1.'.format(vox_number))
		print('')
		sys.exit(1)

	### Make sure numeric inputs are not negative
	if(lag_number<=0):
		print('')
		print('ERROR: the lag number {} must be grater than 0. Exiting with 1.'.format(lag_number))
		print('')
		sys.exit(1)

	if(boot_number<=0):
		print('')
		print('ERROR: the bootstrap number {} must be grater than 0. Exiting with 1.'.format(boot_number))
		print('')
		sys.exit(1)

	if(vox_number<=0):
		print('')
		print('ERROR: the voxel number {} must be grater than 0. Exiting with 1.'.format(vox_number))
		print('')
		sys.exit(1)

	### Make sure numeric inputs are integers
	if(float(int(lag_number)) - lag_number!=0):
		print('')
		print('ERROR: the lag number {} must be an integer. Exiting with 1.'.format(lag_number))
		print('')
		sys.exit(1)

	if(float(int(boot_number)) - boot_number!=0):
		print('')
		print('ERROR: the bootstrap number {} must be an integer. Exiting with 1.'.format(boot_number))
		print('')
		sys.exit(1)

	if(float(int(vox_number)) - vox_number!=0):
		print('')
		print('ERROR: the voxel number {} must be an integer. Exiting with 1.'.format(vox_number))
		print('')
		sys.exit(1)

	### Convert numeric inputs to integers
	lag_number = int(lag_number)
	boot_number = int(boot_number)
	vox_number = int(vox_number)

	### Get image resolution in mm
	input_resmm = input_header.get_zooms()
	input_resmm = np.array(input_resmm)

	### Convert input and ROI to numpy array
	input_data = np.array(input_data)
	roi_data = np.array(roi_data)

	### Binarise the ROI and check for potential errors
	Nvoxel = len(roi_data[roi_data==1])   # Number of voxels in input NIFTI
	if(Nvoxel<2):
		print('ERROR! The ROI mask contains less than 2 voxels and no variogram analysis is possible. Exiting with 1, no output files saved!')
		print('')
		sys.exit(1)
	if(Nvoxel<=vox_number):
		print('WARNING! The ROI contains less voxels than the number of voxels requested for each bootstrap iteration!')
		print('All voxels will be used and only one interation will be performed')
		print('')
		vox_number=Nvoxel
		boot_number=1


	### Get voxel positions in mm and input values as arrays
	print('')	
	print('Allocating data for variogram analysis')
	print('')
	spatialdata = np.zeros((Nvoxel,4))
	voxid = 0
	for ii in range(0,input_size[0]):
		for jj in range(0,input_size[1]):
			for kk in range(0,input_size[2]):

				# Store information regarding a voxel only if it is within the ROI mask
				if(roi_data[ii,jj,kk]==1):
					spatialdata[voxid,0] = input_resmm[0]*ii*1.0
					spatialdata[voxid,1] = input_resmm[1]*jj*1.0
					spatialdata[voxid,2] = input_resmm[2]*kk*1.0
					spatialdata[voxid,3] = input_data[ii,jj,kk]*1.0
					voxid = voxid + 1

	###  Check whether there are at least two different values for the metric whose spatial variability is of interest	
	funique = np.unique(spatialdata[:,3])
	if(funique.size<2):
		print('ERROR! The metric within the ROI contains only one unique value and variogram analysis is not possible. Exiting with 1, no output files saved!')
		print('')
		sys.exit(1)
	

	### Initialise matrices to store fitting output
	ok3d_lags = np.zeros((boot_number,lag_number))
	ok3d_semivar = np.zeros((boot_number,lag_number))
	ok3d_pars =  np.zeros((boot_number,4))
	ok3d_lags[:] = np.nan
	ok3d_semivar[:] = np.nan
	ok3d_pars[:] = np.nan	
	ok3d_pars[:,0] = np.arange(boot_number) + 1

	### Loop over bootstrap calculations
	for bb in range(0,boot_number):

		print('')
		print('     .... variogram fitting {} out of {}. Please wait...'.format(bb+1,boot_number))
		print('')

		## Generate indices of a random sample of voxels if number of voxels to analyse is less than total number of voxels
		if(vox_number<Nvoxel):
			bootidx = np.sort( np.random.choice(Nvoxel,size=vox_number) )
			# Extract the sample
			datasample = spatialdata[bootidx,:]
		## If all voxels need to be analysed, just take them all
		else:
			datasample = np.copy(spatialdata)


		## Check that the sample contains more than one value for the metric of interest
		funiquesample = np.unique(datasample[:,3])
		if(funiquesample.size>=2):
		
			# Yes: the drawn sample contains more than one value for the metric of interest	and variogram analysis is possible
			try:
				# Run variogram analysis
				ok3d = OrdinaryKriging3D(datasample[:,0], datasample[:,1], datasample[:,2], datasample[:,3], variogram_model='exponential', nlags=lag_number)
				
				# Collect output from variogram fitting for this bootstrap iteration
				ok3d_pars[bb,1] = ok3d.variogram_model_parameters[0]
				ok3d_pars[bb,2] = ok3d.variogram_model_parameters[1]
				ok3d_pars[bb,3] = ok3d.variogram_model_parameters[2]
				ok3d_lags[bb,:] = np.array(ok3d.lags)
				ok3d_semivar[bb,:] = np.array(ok3d.semivariance)
			
			# Deal with known errors
			except KeyboardInterrupt:
				print('')
				print('          ERROR: KeyboardInterrupt! Exiting with 1.')
				print('')
				sys.exit(1)	

			except MemoryError:
				print('')
				print('          ERROR: not enough memory! Try to reduce the number of voxels to analyse. Exiting with 1.')
				print('')
				sys.exit(1)

			except np.linalg.linalg.LinAlgError:
				print('')
				print('          WARNING: iteration {} out of {} failed. NaN will be stored in the corresponding entries of the output files.'.format(bb+1,boot_number))
				print('')

			except ValueError:
				print('')
				print('          WARNING: iteration {} out of {} failed. NaN will be stored in the corresponding entries of the output files.'.format(bb+1,boot_number))
				print('')



	
	
	## Save output files with information from all bootstrap calculations
	print('')	
	print('Saving output files')
	print('')
	np.savetxt('{}_semivar.csv'.format(output_base), ok3d_semivar, fmt='%.3e', delimiter=',',comments='')
	np.savetxt('{}_lags.csv'.format(output_base), ok3d_lags, fmt='%.3e', delimiter=',',comments='')
	np.savetxt('{}_fit.csv'.format(output_base), ok3d_pars, fmt='%.3e', delimiter=',', header='bootnr,sill,range,nugget',comments='')
	
	## Done!
	print('')
	print('Processing completed.')
	print('')
	print('')




# Run the module as a script when required
if __name__ == "__main__":

	### Print help and parse arguments
	parser = argparse.ArgumentParser(description='Perform spatial analysis following the VARIO-MS framework (exponential fitting of 3D semivariance and lag from variogram analysis). Authors: Carmen Tur <c.tur@ucl.ac.uk> <carmen.tur.gomez@gmail.com> and Francesco Grussu <f.grussu@ucl.ac.uk> <francegrussu@gmail.com>. # Code released under BSD Two-Clause license # Copyright (c) 2020 University College London. # All rights reserved. ## Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met: ## 1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer. # 2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution. ## THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE. ## The views and conclusions contained in the software and documentation are those of the authors and should not be interpreted as representing official policies, either expressed or implied, of the FreeBSD Project.')
	parser.add_argument('input_nifti', help='input NIFTI whose spatial variability is of interest')
	parser.add_argument('mask_nifti', help='input NIFTI storing a mask with the anatomical region of interest (for instance, the brain)')	
	parser.add_argument('out_basename', help='base name of output files, to which the following suffixes will be added: *_semivar.csv, *_lags.csv, *_fit.csv. *_semivar.csv stores the semivariance across bootstrap iterations (a matrix whose rows are bootstrap iterations and whose columns show values of semivariance for increasing lag). The units of the elements of the matrix are those of the input metric (input_nifti) to the power of two. *_lags.csv stores the lags across bootstrap iterations corresponding to the values of semivariance reported in *_semivar.csv  (it is again a matrix whose rows are bootstrap iterations and whose columns show the lags at which semivariances are calculated). The units of the elements of the matrix are [mm]. *_fit.csv is a CSV file storing the parameters of an exponential fit of the semivariance and lag data. Rows correspond to different bootstrap iterations, while columns correspond to bootstrap index; sill (same units as semivariance), range (in [mm]), nugget (same units as semivariance).')
	parser.add_argument('--lag', metavar='<number>', default='40', help='number of lags for variogram analysis (default: 40)')
	parser.add_argument('--boot', metavar='<number>', default='100', help='number of bootstrap variogram calculations (default: 100)')
	parser.add_argument('--vox', metavar='<number>', default='10000', help='number of voxels subsampled at random to be used for variogram fitting in each bootstrap calculation (default: 10000)')	
	args = parser.parse_args()

	### Get input arguments
	inputfile = args.input_nifti
	maskfile = args.mask_nifti
	outname = args.out_basename	
	lagnr = args.lag
	bootnr = args.boot
	voxnr = args.vox

	### Print a message
	print('')
	print('*************************************************************************************')
	print('                                  VARIO-MS analysis                                  ')
	print('*************************************************************************************')
	print('')
	print('Input NIFTI file: {}'.format(inputfile))
	print('Mask NIFTI gile: {}'.format(maskfile))
	print('Output files: {}_data.csv {}_var.csv {}_fit.csv'.format(outname,outname,outname))
	print('Number of lags: {}'.format(lagnr))
	print('Number of bootstrap calculations: {}'.format(bootnr))
	print('Number of voxels per bootstrap calculation: {}'.format(voxnr))
	print('')

	### Run VARIO-MS analysis
	runVARIOMS(inputfile,maskfile,outname,lagnr,bootnr,voxnr)
	
	### Done: exit with 0
	sys.exit(0)


