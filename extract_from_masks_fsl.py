def extract_from_masks(mask_files,data_files,timeseries=False):
    #given a list of mask files and data files (copes, preprocessed data, etc)
    #returns a data frame of values extracted from those masks
    #the column name is taken based on the mask file in the first row (with file extension stripped)
   
    if len(mask_files) != len(data_files):
        print("mask and data files don't match!")
        
    maskname = mask_files[0].split('/')
    maskname = maskname[len(maskname)-1].strip().split('.')[0]
    
    if timeseries:
        tsdata = []
    else:
        tsdata = pd.Series(data=np.empty(len(mask_files)))
        tsdata[:] = np.NAN

    
    for i,m in enumerate(mask_files):
        ts_out_file = tempfile.NamedTemporaryFile()
        #nipype wrapper for fslmeants
        meants = fsl.ImageMeants()
        meants.inputs.in_file = data_files[i]
        meants.inputs.mask = os.path.abspath(os.path.expanduser(m))
        meants.inputs.out_file = ts_out_file.name
        result = meants.run()
        
        f = open(result.outputs.out_file,'r')
        temp = [line.strip() for line in f.readlines()] #remove the newline character \n
        
        if timeseries:
            tsdata.append(temp)
            print('finished file %s' %data_files[i])
        else:
            tsdata[i] = temp[0]
            
        f.close()
    
    if timeseries:
        output = tsdata
    else:
        output = pd.DataFrame(data=tsdata,columns=[maskname])
    
    return(output)    