from nilearn.masking import apply_mask
import os.path as op

def extract_from_masks(mask_files,data_file):
    #given a data file and a set of mask files, extract timeseries from the data 
    #for each mask. Returns a dictionary with the key corresponding to the mask
    #(derived from the file name) and the value as a numpy array that's TR x voxels
    output = [apply_mask(data_file,m) for m in mask_files] 
    masknames = [op.splitext(op.splitext(op.basename(x))[0])[0] for x in mask_files]
    output = dict(zip(masknames,output))
    return(output) 

def generate_masks(labeled_image_data, therange = None):
    #given an array of image data (e.g., from nib.load(fname).get_data()), find the range of 
    #voxel labels and for each value in the range, add a binary mask with voxels
    #having that value set to TRUE. 
    #Returns the list of (probably) 3d arrays of masks.
    if therange == None:
        therange = range(labeled_image_data.min(), labeled_image_data.max()+1)
    masks = list()
    for i in range(len(therange)):
        masks.append(labeled_image_data == therange[i])
    return(masks)