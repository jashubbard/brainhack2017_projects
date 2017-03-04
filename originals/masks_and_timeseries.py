
# coding: utf-8

# In[2]:

import nipype.interfaces.fsl as fsl
import os
from glob import glob
import tempfile
import numpy as np
import pandas as pd

get_ipython().magic(u'pylab inline')


maskdir = '/Users/canelab/20130226/masks'

os.chdir(maskdir)

if not os.path.exists('individual'):
    os.mkdir('individual')


# In[ ]:




# In[ ]:

#transforms allmasks into native space for each subject
#then extracts timeseries from each mask in the filtered_func data

# dirs = glob('/Users/canelab/20130226/feat_output/firstlevel/000*')
# allsubs=[x.split('/')[6] for x in dirs]

# allmasks= glob('win_charity*')

# for i,subid in enumerate(allsubs):
    #individual subject's feat directory
#     featdir=os.path.join('/Users/canelab/20130226/feat_output/firstlevel/',subid,'r1.feat')
#     regdir = os.path.join(featdir,'reg') #registration directory
#     ind_maskdir=os.path.join(featdir,'masks') #individual mask directory
    
#     print 'starting subject %d of %d:' %(i,len(allsubs))
    
#     if not os.path.exists(ind_maskdir):
#         os.mkdir(ind_maskdir) #make the individual directory if not there
    #for each mask
#     for maskname in allmasks: 
        
#         temp_mat = tempfile.NamedTemporaryFile() #create a temporary file
          
    #nipype wrapper for FSL's applyxfm, warps mask from standard into native space
#         applyxfm = fsl.ApplyXfm()
#         applyxfm.inputs.in_file = maskname
#         applyxfm.inputs.in_matrix_file=os.path.join(regdir,'standard2example_func.mat')
#         applyxfm.inputs.out_file=os.path.join(ind_maskdir,maskname)
#         applyxfm.inputs.out_matrix_file = temp_mat.name
#         applyxfm.inputs.reference = os.path.join(featdir,'example_func.nii.gz')
#         applyxfm.inputs.apply_xfm= True
        
#         result = applyxfm.run()
        
#         print result.outputs
          
          #nipype wrapper for fslmeants
#         meants = fsl.ImageMeants()
#         meants.inputs.in_file = os.path.join(featdir,'filtered_func_data.nii.gz')
#         meants.inputs.mask = os.path.join(ind_maskdir,maskname)
#         meants.inputs.out_file = os.path.join(ind_maskdir,maskname.split('.')[0]+'.txt')
#         result2 = meants.run()
        
#         print result2.outputs


# In[22]:

#general-purpose function for the cell above
def mask_and_extract(featdir,allmasks,img_to_extract):
    
    regdir = os.path.join(featdir,'reg') #registration directory
    ind_maskdir=os.path.join(featdir,'masks') #individual mask directory
    
    
    if not os.path.exists(ind_maskdir):
        os.mkdir(ind_maskdir) #make the individual directory if not there
    
    #for each mask
    for maskname in allmasks: 
        
        mask_out_file = os.path.join(ind_maskdir,maskname);
        
        if not os.path.exists(mask_out_file):
            temp_mat = tempfile.NamedTemporaryFile() #create a temporary file
            #nipype wrapper for FSL's applyxfm, warps mask from standard into native space
            applyxfm = fsl.ApplyXfm()
            applyxfm.inputs.in_file = maskname
            applyxfm.inputs.in_matrix_file=os.path.join(regdir,'standard2example_func.mat')
            applyxfm.inputs.out_file=mask_out_file
            applyxfm.inputs.out_matrix_file = temp_mat.name
            applyxfm.inputs.reference = os.path.join(featdir,'example_func.nii.gz')
            applyxfm.inputs.apply_xfm= True
            
            result = applyxfm.run()
#             print result.outputs

        imgname=img_to_extract.split('.')[0].replace('/','_')
        ts_out_file = img_to_extract.replace('/','_')
        ts_out_file = os.path.join(ind_maskdir,maskname.split('.')[0]+imgname+'.txt')
        
        in_file=featdir+img_to_extract
        
        print "extracting from %s into file %s" %(in_file,ts_out_file)
        
        #nipype wrapper for fslmeants
        meants = fsl.ImageMeants()
        meants.inputs.in_file = in_file
        meants.inputs.mask = os.path.join(ind_maskdir,maskname)
        meants.inputs.out_file = ts_out_file
        result2 = meants.run()
        
#         print result2.outputs

#transforms allmasks into native space for each subject
#then extracts timeseries from each mask in the filtered_func data

dirs = glob('/Users/canelab/20130226/feat_output/firstlevel/000*')
allsubs=[x.split('/')[6] for x in dirs]

os.chdir('/Users/canelab/20130226/masks/harbaugh_mayr_2007')
allmasks= glob('HM*')

for i,subid in enumerate(allsubs):
    #ndividual subject's feat directory
    featdir=os.path.join('/Users/canelab/20130226/feat_output/firstlevel/',subid,'r1.feat')

    
    print 'starting subject %d of %d:' %(i,len(allsubs))
    
    imgnames = ['/stats/cope%d.nii.gz' %num for num in range(1,27)]
    
    for img in imgnames:
            mask_and_extract(featdir,allmasks,img)
#             print featdir,allmasks,img
    print 'complete'



# In[ ]:




# In[15]:

img


# In[19]:

#grabs all the data extracted from the masks and cope images in puts into 1 big dataset
#text files are in the "masks" directory, and are in the format: win_charity_L_nacc_cope13.txt. 
#each contains 1 number. We want 1 row per subject, which a column for each mask and contrast (144 columns)

#get all subject ids
dirs = glob('/Users/canelab/20130226/feat_output/firstlevel/000*')
allsubs=[x.split('/')[6] for x in dirs]

os.chdir('/Users/canelab/20130226/masks/harbaugh_mayr_2007')
allmasks= glob('HM*')

imgnames = ['varcope%d' %num for num in range(1,16)]

alldata= pd.DataFrame(dtype=float)

for i,sub in enumerate(allsubs):
    
    print 'starting subject %d of %d' %(i+1,len(allsubs))
    
    featdir = os.path.join('/Users/canelab/20130226/feat_output/firstlevel/',sub,'r1.feat/')
    ind_maskdir = os.path.join(featdir,'masks')
    
    os.chdir(ind_maskdir)
    
 
    temp = pd.DataFrame(dtype=float)
    
    for img in imgnames:
        tsfiles = glob('HM*_'+img+'.txt')
        
        fnames = [ts.split('.')[0].replace('win_charity_','').replace('stats_','') for ts in tsfiles]
        
        temp2 = pd.DataFrame(data=np.zeros([1,len(fnames)],dtype=float),columns=fnames)
        
        
        for ts in tsfiles:
            f = open(ts)
            data = f.read().split('\n', 1)[0]
            colname = ts.split('.')[0].replace('win_charity_','').replace('stats_','')
            temp2.loc[0,colname] = data
        
        temp = pd.concat((temp,temp2),axis=1)
    
    temp['id'] =sub
   

    alldata = pd.concat((alldata,temp),axis=0)
#   
cols = alldata.columns.tolist()
cols = cols[-1:] + cols[:-1]
alldata = alldata[cols]

alldata.head()


# In[20]:

os.chdir('/Users/canelab/Dropbox/Observed_Altruism/braindata')

alldata.to_csv('HM07_varcope_data.txt',sep='\t',header=True,index=False)


# In[21]:

alldata.columns.tolist()


# In[ ]:




# In[ ]:

allmasks= glob('win_charity*.nii.gz')

alldata = np.empty([0,len(allmasks)+2])

for i,subid in enumerate(allsubs):

    featdir=os.path.join('/Users/canelab/20130226/feat_output/firstlevel/',subid,'r1.feat')
    regdir = os.path.join(featdir,'reg')
    ind_maskdir=os.path.join(featdir,'masks')
    
    os.chdir(ind_maskdir)
    
    temp = pd.read_csv(allmasks[0].split('.')[0]+'.txt',header=None)
    nrows=len(temp)
    
    subdata=np.zeros((nrows,len(allmasks)+2))
    
    evs = np.loadtxt(featdir+'/ev_tr.txt')
    
    subdata[:,0] = np.repeat(subid,nrows)
    subdata[:,1] = evs[:,0]
    
    for j,m in enumerate(allmasks):
        maskname=m.split('.')[0]
        
        ts = np.loadtxt(maskname+'.txt')
        
        subdata[:,j+2] = ts
    
    alldata=np.concatenate((alldata,subdata))
    print alldata.shape


# In[ ]:

maskheadernames = [s.split('.')[0] for s in allmasks]


finaldata = pd.DataFrame(alldata,columns=['id','condition']+maskheadernames)

#add leading zeros to id number
temp = finaldata['id'].astype(int).values
finaldata['id'] = ['%07d' %i for i in temp]



os.chdir(maskdir)
finaldata.to_csv('all_masked_timeseries.txt',sep='\t',header=True,index=False)



# In[ ]:

def ev_to_tr(fsfdir):
#     print "Starting subject %s:\n" %s
    
    import numpy as np
    from glob import glob
    from mvpa2.misc.fsl.base import read_fsl_design
    from mvpa2.misc.fsl.base import FslGLMDesign
    import os 
    
    
    if os.path.exists(os.path.join(fsfdir,'design.mat')):
    
        #fsfdir=os.path.join(rootdir,s,featdir)
        #print fsfdir
        fsffile=''.join([fsfdir,'design.fsf'])
        desmatfile=''.join([fsfdir,'design.mat'])
    
        design=read_fsl_design(fsffile)
        
        desmat=FslGLMDesign(desmatfile)
        
        
        ntp=desmat.mat.shape[0]
    
        TR=design['fmri(tr)']
        
        #print "EVs: %d" %nevs
        #print "timepoints: %d" %ntp
        #print "TR: %d" %TR
        
        evdir=os.path.join(fsfdir,"custom_timing_files")
        os.chdir(evdir)
        evfiles=glob('ev*.txt')
        nevs = len(evfiles)
        #print evfiles
        
        #where we'll code the data
        alldata=np.zeros([ntp, len(evfiles)+1],dtype='float')
        
        allevdata=np.zeros([1,3],dtype='float')
        
        count=0
        for e in evfiles:
            #print count+1
            evdata=np.loadtxt(open(e,"rb"),delimiter="\t",skiprows=0,dtype='float')
            trdata=np.int_(np.ceil(evdata/TR)) 
          
            temp = np.zeros([len(evdata),3],dtype="float")
            temp[:,0] = np.repeat(count+1,len(temp))
            temp[:,1:4] = evdata[:,0:2]
            
            if count==0:
                allevdata=temp
            else:
                allevdata = np.vstack((allevdata,temp))
            
            ind= np.lexsort((allevdata[:,0],allevdata[:,1]))
            allevdata = allevdata[ind];
            
            durs=trdata[:,1]
            
            #durs=durs[durs>=0]
            
            #print evdata[:,1]
            onsets=trdata[:,0]
            endtimes=onsets+durs
            endtimes[endtimes>ntp]=ntp
            intervals=np.transpose(np.vstack((onsets,endtimes)))
    
            intervals=intervals[intervals[:,0]<ntp,:]
            intervals=intervals[intervals[:,0]>0,:]
            
            for i in range(intervals.shape[0]):
                trs=np.arange(intervals[i,0],intervals[i,1]+1)
                
                trs=trs[trs<ntp]
                #print trs
                alldata[trs,0]=count+1
                alldata[trs,count+1]=1
            
            
            
            count +=1
        evtrfile=os.path.join(fsfdir,'ev_tr.txt')
        allevfile=os.path.join(fsfdir,'allevs.txt')
        print "saving %s" %evtrfile
        np.savetxt(evtrfile,alldata,delimiter='\t',fmt='%d')
        print "saving %s" %allevfile
        np.savetxt(allevfile,allevdata,delimiter='\t',fmt=['%d','%8.3f','%8.3f'])
       
        
    
        


        
    else:
        print "skipping directory %s. No design.mat file!" %fsfdir 


# In[ ]:




# In[ ]:

for subdir in dirs:
    ev_to_tr(subdir+'/r1.feat/')


# In[ ]:

#takes a list of onsets and durations (in seconds) and expands into 
#a tr-by-tr list, fills in "val" at appropriate rows

def expand_to_tr(onsets,durs,TR,ntp,val):
    
    alltrs = np.zeros([ntp,1],dtype=int)
    
    onsets = np.ceil(onsets/TR).astype(int)
    durs = np.ceil(durs/TR).astype(int)
    endtimes=onsets+durs
    
    
    endtimes[endtimes>ntp]=ntp
    intervals=np.transpose(np.vstack((onsets,endtimes))-1)

    intervals=intervals[intervals[:,0]<ntp,:]
    intervals=intervals[intervals[:,0]>0,:]

    
    for i in range(intervals.shape[0]):
        trs=np.arange(intervals[i,0],intervals[i,1]+1)
        trs=trs[trs<ntp]
        alltrs[trs]=val
        
    return alltrs
    
    


# In[ ]:

test = np.array(ev_subset.onset)
test2 = np.ceil(test/TR).astype(int)


print test,test2


# In[ ]:

#look at all evs and convert to tr-level data (and adjust time periods)
#then save in a dataset with subid, evnumber, then timeseries from each mask

from mvpa2.misc.fsl.base import FslGLMDesign #for reading in the design.mat file to get the number of timepoints

dirs = glob('/Users/canelab/20130226/feat_output/firstlevel/000*')
allsubs=[x.split('/')[6] for x in dirs] #get all the subject ids from the feat directory names


maskdir = '/Users/canelab/20130226/masks'
os.chdir(maskdir)
# alldirs = [d+'/r1.feat/' for d in dirs];

#the ev numbers that we want to extract from
evnums = [5,6,7,8]
TR = 2.

onset_adjust=4.
window_dur=6.

allmasks= glob('win_charity*.nii.gz')
alldata = np.empty([0,len(allmasks)+2])

for i,subid in enumerate(allsubs):
    
    featdir=os.path.join('/Users/canelab/20130226/feat_output/firstlevel/',subid,'r1.feat/')
    desmat=FslGLMDesign(featdir+'design.mat')
    ntp=desmat.mat.shape[0]
    
#     regdir = os.path.join(featdir,'reg')
    ind_maskdir=os.path.join(featdir,'masks') #mask directory in individual feat dir
    os.chdir(ind_maskdir)
    
#     temp = pd.read_csv(allmasks[0].split('.')[0]+'.txt',header=None)
#     nrows=len(temp)
    
    subdata=np.zeros((ntp,len(allmasks)+2)) #this is where the each subject's data goes
    
    #read in all evs 
    allevs=pd.read_csv(featdir+'allevs.txt',sep='\t',header=None,names=['ev','onset','duration'])
    
    #after converting to trs, you will save it here
    #create 1 column for each ev
    alltrdata = np.zeros((ntp,len(evnums)),dtype=int)
    
    #loop through each ev specified above
    for i,evnum in enumerate(evnums):    
        
        #grab a subset from the larger dataset
        ev_subset = allevs.loc[allevs.ev==evnum,:]

        #adjust timing and duration to reflect the window you want around stimulus onset
        ev_subset.onset = ev_subset.onset+onset_adjust
        ev_subset.duration[:] = window_dur
        
        #expand them out to trs, then save as a column in our dataset
        alltrdata[:,i] = expand_to_tr(ev_subset.onset,ev_subset.duration,TR,ntp,evnum)[:,0]
    
    #fill in subdata with subject name, then the tr information (collapsed across columns)
    subdata[:,0] = np.repeat(subid,ntp)
    subdata[:,1] = np.sum(alltrdata,1)
    
    
    #loop through each mask 
    for j,m in enumerate(allmasks):
        maskname=m.split('.')[0]
        ts = np.loadtxt(maskname+'.txt') #load the timeseries
        subdata[:,j+2] = ts #store in subject data in appropriate column (starting with 3rd column)
    
    #add the subject data to master dataset
    alldata=np.concatenate((alldata,subdata))
    print 'subject %s complete' %subid

    
  



# In[ ]:

#puts all the data in a pandas dataframe and exports to csv
colnames=['id','ev']+[fname.split('.')[0] for fname in allmasks]

finaldata = pd.DataFrame(alldata,columns=colnames)



#add leading zeros to id number
temp = finaldata['id'].astype(int).values
finaldata['id'] = ['%07d' %i for i in temp]

os.chdir('/Users/canelab/Dropbox/Observed_Altruism/braindata')
finaldata.to_csv('allmaskdata_adjusted_Sc2.txt',sep='\t',index=False)


# In[ ]:

#start here after saving everything


maskdir = '/Users/canelab/20130226/masks'
os.chdir(maskdir)

df = pd.read_csv('allmaskdata.txt',sep='\t',dtype={'id':str})

df = df.loc[df.ev>0,:]

df.head()


            


# In[ ]:

grouped = df.groupby(['id','ev'])
aggsub = grouped.aggregate(np.mean)


# In[ ]:




# In[ ]:




# In[ ]:



