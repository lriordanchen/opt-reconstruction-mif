"""
This notebook reads in hdf5 files and constructs tomograms 
using the TomoPy package 
( https://tomopy.readthedocs.io/en/latest/about.html ).

You may need to install some packages to get this to work
( https://jakevdp.github.io/blog/2017/12/05/installing-python-packages-from-jupyter/ ). 
"""
import sys
import tomopy
import matplotlib.pyplot as plt
import h5py
import imageio
import numpy as np
import pickle
import datetime
from skimage import transform as transf
from skimage import img_as_float64
import os
from tifffile import imsave

def opt_recon(dirname, tilt_correction=False, mlem_option=False):
	# Find all files in input 
	dirname_import = os.path.join(dirname,r'input')
	print(dirname)
	file_set = []
	for path, dirs, files in os.walk(dirname_import):
		for f in files:
			file_set.append(os.path.join(path,f))


	file_set = [s for s in file_set if '.h5' in s]
	
	# Grab names of each file, represents different channels
	filter_start = '\\input\\'
	filter_end = '.h5'
	file_names = [s[s.find(filter_start)+len(filter_start):s.find(filter_end)] for s in file_set]
	
	# Check for background and dark-field files
	flat_bool=False
	dark_bool=False
	if any('bkgd' in s  for s in file_names):
		flat_bool=True
	if any('dark' in s  for s in file_names):
		flat_bool=True
	
	# Put transmission file first
	if ('trans' in s for s in file_names):
		trans_bool=True
		fnames = [dirname_import+r'\trans.h5']
		file_names_add = [s for s in file_names if not 'trans' in s]
	else:
		trans_bool=False
		fnames = []
	# Add the rest of the files to fnames for reading in
	for ix in range(len(file_names_add)):
		fnames.append(os.path.join(dirname_import, file_names_add[ix])+'.h5')
	print(fnames)


	# Read hdf5 file structure. 
	# This extracts the structure in the hdf5 file (names of keys, 
	# number of images, size of images, data type) and displays it below.		 
	dset_array = np.array(range(len(fnames)), dtype='|S20')
	for ix in range(len(fnames)):
		filename = fnames[ix]
		f = h5py.File(filename, 'r')


		def traverse_datasets(hdf_file):

			def h5py_dataset_iterator(g, prefix=''):
				for key in g.keys():
					item = g[key]
					path = '{0}/{1}'.format(prefix, key)
					if isinstance(item, h5py.Dataset): # test for dataset
						yield (path, item)
					elif isinstance(item, h5py.Group): # test for group (go down)
						for ix in h5py_dataset_iterator(item, path):
							yield ix

			with h5py.File(hdf_file, 'r') as f:
				for path, _ in h5py_dataset_iterator(f):
					yield path
				
		with h5py.File(filename, 'r') as f:
			for dset in traverse_datasets(filename):
				print('Path:', dset)
				dset_array[ix]=dset
				print('Shape:', f[dset].shape)
				shape_var = f[dset].shape
				len_var= sum(len(x) for x in f[dset])
				print('Data type:', f[dset].dtype)		

	# Read and save the data from the hdf5 file in variable *proj*
	proj = np.zeros(shape=(len(fnames), shape_var[0], shape_var[1], shape_var[2]))

	for ix in range(len(fnames)):
		print('Reading out channel: ',ix+1)
		filename = fnames[ix]
		f = h5py.File(filename, 'r')
		dataset=f[dset_array[ix]]
		proj[ix] = np.array(dataset[:,:,:])
	print(proj.shape)
	num_channels, num_images, image_height, image_width = proj.shape
	print('Import hdf5 complete.')		

	# Flat-field and dark-field corrections	

	if flat_bool==True:
		filename_flat = dirname+r'\input\trans_bkgd.tif'
		flat = plt.imread(filename_flat)

	if dark_bool==True:
		filename_dark = r'\input\trans_dark.tif'
		dark = plt.imread(filename_dark)

	if flat_bool==True and dark_bool==False:
		#We don't often take dark-field images...
		dark_value = 0
		dark = np.full((proj.shape[2],proj.shape[3]), dark_value)

	if flat_bool==True:
		proj[0] = tomopy.normalize(proj[0], flat, dark)
		
	#The following assumes 360 degrees rotation divided evenly between all images.
	theta = np.linspace(0,np.pi*2,num_images)
	image_180 = int(np.floor(num_images/2)) #assuming 360 degrees of rotation, this is image 180 degrees apart from first image
	try:
		rot_file = open(dirname+r'\input\rotation.pkl','rb')
		rot_center = pickle.load(rot_file)
		rot_file.close()
		print('Rotation file found')
	except FileNotFoundError:
		rot_center = np.zeros(num_channels)
		for ch in range(num_channels):
			rot_center[ch] = tomopy.find_center_pc(proj[ch][0],proj[ch][image_180])
			print(rot_center[ch])
		#Pickle for future reference
		rot_file = open(dirname+r'\input\rotation.pkl','wb')
		pickle.dump(rot_center, rot_file)
		rot_file.close()
	print('Rotation complete.')



	if tilt_correction==True:
		alignment_channel = 1
		#First check if alignment shifts have previously been calculated
		try:
			tilt_file = open(dirname+r'\input\tilt_shift.pkl', 'rb')
			sy, sx = pickle.load(tilt_file)
			tilt_file.close()
			print('Tilt file found')

		#Otherwise run alignment module
		except FileNotFoundError:
			print('Tilt file not found. Using TomoPy to calculate alignment.')
			
			# This can take a long time (~5-10 min/iteration).  Often times, ~10 iterations is sufficient.  Ideally err<1.

			align_params = {'algorithm':'mlem', 'iters':20}

			start = datetime.datetime.now()
			print("Aligning channel: ", alignment_channel)
			proj_align_channel, sy , sx, conv_channel = tomopy.prep.alignment.align_joint(proj[alignment_channel], theta, 
																						  center=rot_center[alignment_channel], 
																						  **align_params)
			end = datetime.datetime.now()
			print( int((end - start).total_seconds()/60), 'minutes' )

			#Pickle for future reference
			tilt_file = open(dirname+r'\input\tilt_shift.pkl','wb')
			pickle.dump([sy,sx], tilt_file)
			tilt_file.close()
		print("Tilt calculation complete.")
	else:
			print('No alignment correction')

	if tilt_correction==True:
		align_proj = np.zeros(shape=proj.shape)
		for ch in range(num_channels):
			print('Aligning channel: ', ch+1)
			channel_proj = proj[ch]
			for ix in range(num_images):
				tform = transf.SimilarityTransform(translation=(sy[ix], sx[ix]))
				channel_proj[ix] = transf.warp(channel_proj[ix], tform)
			align_proj[ch] = channel_proj
		print('Tilt application complete.')
	else:
		print('No alignment correction')
		align_proj=proj
		
		
	# If applied tilt, save new projection -- again, if tilt_projection = False, this does not change anything.
	proj=align_proj

	# Calculate $$ inv(proj) $$
	if trans_bool==True:
		proj_inv = np.zeros(shape = proj[0].shape)
		proj_inv = np.max(proj[0])-proj[0]
		#proj_inv = -np.log(proj)
		align_proj[0]=proj_inv

	
	# Reconstructino using gridrec
	params = {'algorithm':'gridrec', 'filter_name':'butterworth' }

	recon = np.zeros(shape=(num_channels, image_height, image_width, image_width))
	start = datetime.datetime.now()
	for ix in range(num_channels):
		print('Reconstructing channel: ', ix+1)
		recon_channel = tomopy.recon(align_proj[ix], theta, center=rot_center[ix], 
									 **params
									  )
		recon[ix] = recon_channel
	recon[recon<0] = 0
	end = datetime.datetime.now()
	print( int((end - start).total_seconds()/60), 'minutes' )
	print('gridrec reconstruction complete.')

	# Mask with circle
	for ix in range(len(recon)):
		recon[ix] = tomopy.circ_mask(recon[ix], axis=0, ratio=0.95)

	# Save data as tiff and hdf5
	#Set folder and names of outputs, must be same number as number of inputs, ideally use same labels
	dirOutput = dirname + r'\output_gridrec'
	output_filenames = [s[s.find(filter_start)+len(filter_start):s.find(filter_end)] for s in fnames]
	print(output_filenames)

	file_output = []
	for ix in range(len(output_filenames)):
		file_output.append(os.path.join(dirOutput, output_filenames[ix]))

	try:
		# Creates directory
		os.mkdir(dirOutput)
		print("Directory gridrec created ") 
	except FileExistsError:
		print("Directory gridrec already exists")
		
	for ix in range(len(file_output)):
		print('Writing channel: ', ix+1)
		file_output_h5 = file_output[ix]+'.h5'
		archive = h5py.File(file_output_h5, 'w')
		archive['recon'] = recon[ix]
		archive.close()
		
		imsave(file_output[ix]+'.tif', np.array(recon[ix], dtype='float32'))
		
	print('Writing of gridrec output complete.')
	
	
	
	if mlem_option==True:
		# Reconstruction using mlem
		params = {'algorithm':'mlem', 'num_iter':20 }
		recon = np.zeros(shape=(num_channels, image_height, image_width, image_width))
		start = datetime.datetime.now()
		for ix in range(num_channels):
			print('Reconstructing channel: ', ix+1)
			recon_channel = tomopy.recon(align_proj[ix], theta, center=rot_center[ix], 
										 **params
										  )
			recon[ix] = recon_channel
		recon[recon<0] = 0
		end = datetime.datetime.now()
		print( int((end - start).total_seconds()/60), 'minutes' )
		print('mlem reconstruction complete.')

		# Mask with circle
		for ix in range(len(recon)):
			recon[ix] = tomopy.circ_mask(recon[ix], axis=0, ratio=0.95)

		# Save data as tiff and hdf5
		#Set folder and names of outputs, must be same number as number of inputs, ideally use same labels
		dirOutput = dirname + r'\output_mlem'
		output_filenames = [s[s.find(filter_start)+len(filter_start):s.find(filter_end)] for s in fnames]
		print(output_filenames)

		file_output = []
		for ix in range(len(output_filenames)):
			file_output.append(os.path.join(dirOutput, output_filenames[ix]))

		try:
			# Creates directory
			os.mkdir(dirOutput)
			print("Directory mlem created ") 
		except FileExistsError:
			print("Directory mlem already exists")
			
		for ix in range(len(file_output)):
			print('Writing channel: ', ix+1)
			file_output_h5 = file_output[ix]+'.h5'
			archive = h5py.File(file_output_h5, 'w')
			archive['recon'] = recon[ix]
			archive.close()
			
			imsave(file_output[ix]+'.tif', np.array(recon[ix], dtype='float32'))
			
		print('Writing of mlem output complete.')
	
	


if __name__ == "__main__":
	input = sys.argv[1:]
	print(input)
	tilt_corr_bool=False
	mlem_bool=False
	
	
	tilt_corr_bool_string = [s for s in input if 'tilt_correction' in s]
	input = [s for s in input if not 'tilt_correction' in s]
	tilt_corr_bool=True if 'true' in str(tilt_corr_bool_string) else False
	
	mlem_bool_string = [s for s in input if 'mlem' in s]
	input = [s for s in input if not 'mlem' in s]
	mlem_bool=True if 'true' in str(mlem_bool_string) else False
	
	program_start = datetime.datetime.now()
	for ix in range(len(input)):
		opt_recon(input[ix], tilt_correction = tilt_corr_bool, mlem_option = mlem_bool)
	program_end = datetime.datetime.now()
	print( int((program_end - program_start).total_seconds()/60), 'minutes' )