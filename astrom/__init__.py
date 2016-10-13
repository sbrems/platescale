import os
from astropy.io import fits
import numpy as np


#home_dir='/home/sbrems/Documents/NACO/astrometry/'


def make_dirs(save_directories,delete_old=False):
    #check if folders exist and create
    #save_directories = [dir_interm]
    for director in save_directories:
        if not os.path.exists(director):
            os.makedirs(director)
    #check folders are empty and delete old files
    if delete_old:
        print 'Deleting old output files!' 
        for director in save_directories:
            for the_file in os.listdir(director):
                file_path = os.path.join(director, the_file)
                if os.path.isfile(file_path):
                    os.remove(file_path)
                #elif os.path.isdir(file_path): shutil.rmtree(file_path) #del also subfolders


def read_fits(directory):
    '''This routine reads all fits files data into a big data cube and all header files
    into a big header cube. The order is the same and is alphabetically.Filenames and headers
    are multiple if it was an imagecube. So all have the same length. So it returns:
    (fits)filenames,datacube,headercube'''
    #to avoid buffer overflows we need the number images first, which is not the same as
    #filenumber, as some images are cubes and some are not
    n_images = 0
    form = []
    for file in sorted(os.listdir(directory)):
        if file.endswith('.fits'):
            form = fits.getdata(directory+file).shape
            if len(form) == 3 : #image cube
                n_images += form[0]
            elif len(form) == 2: #one image
                n_images += 1
            else:
                raise ValueError('Fits file has unknown format!')
    print('Found ',n_images,' frames in ',directory)
    
    #now make the array
    filenames = []
    headers = []
    all_data = np.full((n_images,form[-2],form[-1]),np.nan,dtype=np.float32) #float16 for memory
    n = 0
    for file in sorted(os.listdir(directory)):
        if file.endswith('.fits'):
            data,header = fits.getdata(directory+file,header=True)
            if len(data.shape) == 3 :#image cube
                all_data[n:n+data.shape[0],:,:] = data
                headers.extend(header for ii in xrange(data.shape[0]))
                filenames.extend(file for ii in xrange(data.shape[0]))
                n += data.shape[0]
            elif len(data.shape) == 2: #one image
                all_data[n,:,:] = data
                headers.append(header)
                filenames.append(file)
                
            else:
                raise ValueError('Fits file has unknown format!')

    return filenames, all_data,headers
