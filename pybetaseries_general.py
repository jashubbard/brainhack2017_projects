#!/usr/bin/env python
"""pybetaseries: a module for computing beta-series regression on fMRI data

Includes:
pybetaseries: main function
spm_hrf: helper function to generate double-gamma HRF
"""

from __future__ import division
# from mvpa.misc.fsl.base import *
from mvpa2.misc.fsl.base import *
from mvpa2.datasets.mri import fmri_dataset
from nilearn.masking import apply_mask
import numpy as np
import nibabel as nib
import os.path as op
import nilearn
import scipy.stats
from scipy.ndimage import convolve1d
from scipy.sparse import spdiags
from scipy.linalg import toeplitz
from mvpa2.datasets.mri import *
import os, sys
from copy import copy
import argparse
import pandas as pd

class pybetaseries():
    
    def __init__(self,img = None, mask=None, design= None, fsl=False, fsfdir=None,outdir = None, tempderiv = None, motpars=None, time_res = None,TR=None,HPF=None):
        #Sets up dataset wide variables
        
        self.tempderiv=tempderiv
        self.motpars=motpars
        self.time_res=time_res
        self.fsl = fsl
        

        if fsl and not os.path.exists(fsfdir):
            print('ERROR: %s does not exist!'%fsfdir)
       
        if fsl:     
            if not fsfdir.endswith('/'):
                fsfdir=''.join([fsfdir,'/'])
            
            self.fsfdir=fsfdir

            fsffile=''.join([self.fsfdir,'design.fsf'])
            desmatfile=''.join([self.fsfdir,'design.mat'])
    
            design=read_fsl_design(fsffile)
    
            self.desmat=FslGLMDesign(desmatfile)
        
            self.nevs=self.desmat.mat.shape[1]
            self.ntp=self.desmat.mat.shape[0]
        
            self.TR=round(design['fmri(tr)'],2)
    
            self.hrf=spm_hrf(self.time_res)

            if not os.path.exists(fsfdir+'betaseries'):
                os.mkdir(fsfdir+'betaseries')

            maskimg=''.join([fsfdir,'mask.nii.gz'])
            self.raw_data=fmri_dataset(fsfdir+'filtered_func_data.nii.gz',mask=maskimg)
            
            voxmeans = np.mean(self.raw_data.samples,axis=0)
            self.data=self.raw_data.samples-voxmeans
            self.nvox=self.raw_data.nfeatures
        	
#             imgdata = nib.load(fsfdir+'filtered_func_data.nii.gz')
            self.mask = nib.load(maskimg)
#             self.raw_data = imgdata
            self.masked_data = []
            cutoff=design['fmri(paradigm_hp)']/self.TR
            self.outdir = op.expanduser(self.fsfdir)

        else:

            self.outdir = op.expanduser(outdir)
            self.process_design(design)

            imgdata = nib.load(img)
            self.mask  = nib.load(mask)
            self.raw_data = imgdata
            self.masked_data = apply_mask(self.raw_data,self.mask)
            self.ntp = imgdata.shape[3]
            self.nevs = self.design.conditions.unique()
            self.TR=round(TR,2)
            self.hrf=spm_hrf(self.time_res)
            cutoff= HPF/self.TR
            
            voxmeans = np.mean(self.masked_data,axis=0)
            self.data=self.masked_data-voxmeans
            self.nvox=self.raw_data.shape[1]
    
    
        self.time_up=np.arange(0,self.TR*self.ntp+self.time_res, self.time_res);
        
        self.max_evtime=self.TR*self.ntp - 2;
        self.n_up=len(self.time_up)
        self.F=get_smoothing_kernel(cutoff, self.ntp)
            
    def LSS(self,whichevs=None,numrealev=None):
        method='LSS'
        # print("Calculating LSS: Ev %s" %(whichevs[0]))
        print('Calculating LSS')
        
        if self.fsl:
            nuisance = otherevs(whichevs,numrealev,self.tempderiv,self.motpars)    
            ons=FslEV3(self.fsfdir+'custom_timing_files/ev%d.txt'%int(whichevs[0]))
            dm_nuisanceevs = self.desmat.mat[:, nuisance]
            
        else:
            ons = self.design.loc[:,['onsets','durations','amplitudes']]
            dm_nuisanceevs = None

            if self.motpars:
                dm_nuisanceevs = self.motpars
        

        ntrials=len(ons.onsets)
        beta_maker=np.zeros((ntrials,self.ntp))
        dm_trials=np.zeros((self.ntp,ntrials))
        
        for t in range(ntrials):
            if ons.onsets[t] > self.max_evtime:
                continue
            # build model for each trial
            dm_trial=np.zeros(self.n_up)
            window_ons = [np.where(self.time_up==x)[0][0]
                      for x in self.time_up
                      if ons.onsets[t] <= x < ons.onsets[t] + ons.durations[t]]
            
            dm_trial[window_ons]=1
            dm_trial=np.convolve(dm_trial,self.hrf)[0:int(self.ntp/self.time_res*self.TR):int(self.TR/self.time_res)]
            dm_trials[:,t]=dm_trial
            print('Finished trial %d of %d' %(t+1,ntrials))

        dm_full=np.dot(self.F,dm_trials)
        dm_full=dm_full - np.kron(np.ones((dm_full.shape[0],dm_full.shape[1])),np.mean(dm_full,0))[0:dm_full.shape[0],0:dm_full.shape[1]]
        
        
            
        for p in range(len(dm_full[1,:])):
            target=dm_full[:,p]
            dmsums=np.sum(dm_full,1)-dm_full[:,p]

            if dm_nuisanceevs is not None:
                des_loop=np.hstack((target[:,np.newaxis],dmsums[:,np.newaxis],dm_nuisanceevs))
            else:
                des_loop=np.hstack((target[:,np.newaxis],dmsums[:,np.newaxis]))

            beta_maker_loop=np.linalg.pinv(des_loop)
            beta_maker[p,:]=beta_maker_loop[0,:]
                
                # des_loop=np.hstack((target[:,np.newaxis],dmsums[:,np.newaxis],dm_nuisanceevs))
        # this uses Jeanette's trick of extracting the beta-forming vector for each
        # trial and putting them together, which allows estimation for all trials
        # at once
        
        glm_res_full=np.dot(beta_maker,self.data)
        # ni=map2nifti(self.raw_data,data=glm_res_full)
        
        # if self.fsl:
        #     ni.to_filename(self.outdir+'betaseries/ev%d_%s.nii.gz'%(int(whichevs[0]),method))
        # else:
        #     ni.to_filename(op.join(self.outdir,self.design.filename[i]))

        if self.fsl:
        	ni=map2nifti(self.raw_data,data=glm_res_full)
        	ni.to_filename(self.outdir+'betaseries/ev%d_%s.nii.gz'%(int(whichevs[0]),method))
        else:
            for t in range(glm_res_full.shape[0]):
                template = nilearn.image.new_img_like(self.mask,np.zeros(self.mask.shape))
                template.get_data()[self.mask.get_data()==1] = glm_res_full[t,:]
                template.to_filename(op.join(self.outdir,self.design.filename[t]+'.nii.gz'))

            
    def LSA(self,whichevs,numrealev):
        method='LSA'
        print("Calculating LSA")
        
        dm_nuisanceevs = None
       
        
        if self.fsl:
            nuisance = otherevs(whichevs,numrealev,self.tempderiv,self.motpars)
            dm_nuisanceevs = self.desmat.mat[:, nuisance]
            all_onsets=[]
            all_durations=[]
            all_conds=[]  # condition marker

            for e in range(len(whichevs)):
                ev=whichevs[e]
                ons=FslEV3(self.fsfdir+'custom_timing_files/ev%d.txt'%int(ev))       
                all_onsets=np.hstack((all_onsets,ons.onsets))
                all_durations=np.hstack((all_durations,ons.durations))
                all_conds=np.hstack((all_conds,np.ones(len(ons.onsets))*(ev)))
        else:
            all_onsets = self.design.onsets
            all_durations = self.design.durations
            all_conds = self.design.conditions


        ntrials=len(all_onsets)
        glm_res_full=np.zeros((self.nvox,ntrials))
        dm_trials=np.zeros((self.ntp,ntrials))
        dm_full=[]
        for t in range(ntrials):
            if all_onsets[t] > self.max_evtime:
                continue
            dm_trial=np.zeros(self.n_up)
            window_ons = [np.where(self.time_up==x)[0][0]
                      for x in self.time_up
                      if all_onsets[t] <= x < all_onsets[t] + all_durations[t]]
            
            dm_trial[window_ons]=1
            dm_trial=np.convolve(dm_trial,self.hrf)[0:int(self.ntp/self.time_res*self.TR):int(self.TR/self.time_res)]
            dm_trials[:,t]=dm_trial
            
        dm_full=np.dot(self.F,dm_trials)
        dm_full=dm_full - np.kron(np.ones((dm_full.shape[0],dm_full.shape[1])),np.mean(dm_full,0))[0:dm_full.shape[0],0:dm_full.shape[1]]
        
  
        if dm_nuisanceevs is not None:        
            dm_full=np.hstack((dm_full,dm_nuisanceevs))
        

        glm_res_full=np.dot(np.linalg.pinv(dm_full),self.data)
        glm_res_full=glm_res_full[0:ntrials,:]
    
        # for e in whichevs:
        #     ni=map2nifti(self.raw_data,data=glm_res_full[np.where(all_conds==(e))[0],:])
        #     ni.to_filename(self.fsfdir+'betaseries/ev%d_%s.nii.gz'%(int(e),method))

        if self.fsl:
            for e in whichevs:
                ni=map2nifti(self.raw_data,data=glm_res_full[np.where(all_conds==(e))[0],:])
                ni.to_filename(self.fsfdir+'betaseries/ev%d_%s.nii.gz'%(int(e),method))
        else:
            for t in range(glm_res_full.shape[0]):
                template = nilearn.image.new_img_like(self.mask,np.zeros(self.mask.shape))
                template.get_data()[self.mask.get_data()==1] = glm_res_full[t,:]
                template.to_filename(op.join(self.outdir,self.design.filename[t]+'.nii.gz'))

    def process_design(self,dfile):

        df = pd.DataFrame.from_csv(dfile,index_col=None)
        df = df.sort_values('onsets')

        if 'amplitudes' not in df:
            df.loc[:,'amplitudes'] = np.ones([len(df.onsets),1])

        
        conds = df['conditions'].unique()


        df.loc[:,'trial'] = np.arange(1,len(df.onsets)+1)
        fnames = ['%03d_%s' % (row['trial'],row['conditions']) for (i,row) in df.loc[:,['trial','conditions']].iterrows()]
        df.loc[:,'filename'] = fnames
        df = df[['trial','conditions','onsets','durations','amplitudes','filename']]
        
        #write it 
        df.to_csv(op.join(self.outdir,'design.csv'),index=False,header=True)
        self.design = df

        return(df)

        
def get_smoothing_kernel(cutoff, ntp):
    sigN2 = (cutoff/(np.sqrt(2.0)))**2.0
    K = toeplitz(1
                 /np.sqrt(2.0*np.pi*sigN2)
                 *np.exp((-1*np.array(range(ntp))**2.0/(2*sigN2))))
    K = spdiags(1./np.sum(K.T, 0).T, 0, ntp, ntp)*K
    H = np.zeros((ntp, ntp)) # Smoothing matrix, s.t. H*y is smooth line
    X = np.hstack((np.ones((ntp, 1)), np.arange(1, ntp+1).T[:, np.newaxis]))
    for  k in range(ntp):
        W = np.diag(K[k, :])
        Hat = np.dot(np.dot(X, np.linalg.pinv(np.dot(W, X))), W)
        H[k, :] = Hat[k, :]

    F = np.eye(ntp) - H
    return F

def otherevs(whichevs,numrealev,tempderiv,motpars):
        #sets up the onsets and nuisance EVs for given target EV 
        
        if tempderiv:
            nuisance=range(0,2*numrealev)
            popevs=[(ev-1)*2 for ev in whichevs]
            nuisance=[i for i in nuisance if i not in popevs]

        
            if motpars:
                nuisance.extend(range(2*numrealev,(6*2+2*numrealev)))
        
        
        else:
            nuisance=range(0,numrealev)
            popevs=[(ev-1) for ev in whichevs]
            nuisance=[i for i in nuisance if i not in popevs]
        
            if motpars:
                nuisance.extend(range(numrealev,6+numrealev))
        
        return nuisance



def spm_hrf(TR,p=[6,16,1,1,6,0,32]):
    """ An implementation of spm_hrf.m from the SPM distribution

    Arguments:

    Required:
    TR: repetition time at which to generate the HRF (in seconds)

    Optional:
    p: list with parameters of the two gamma functions:
                                                         defaults
                                                        (seconds)
       p[0] - delay of response (relative to onset)         6
       p[1] - delay of undershoot (relative to onset)      16
       p[2] - dispersion of response                        1
       p[3] - dispersion of undershoot                      1
       p[4] - ratio of response to undershoot               6
       p[5] - onset (seconds)                               0
       p[6] - length of kernel (seconds)                   32

    """

    p=[float(x) for x in p]

    fMRI_T = 16.0

    TR=float(TR)
    dt  = TR/fMRI_T
    u   = np.arange(p[6]/dt + 1) - p[5]/dt
    hrf=scipy.stats.gamma.pdf(u,p[0]/p[2],scale=1.0/(dt/p[2])) - scipy.stats.gamma.pdf(u,p[1]/p[3],scale=1.0/(dt/p[3]))/p[4]
    good_pts=np.array(range(np.int(p[6]/TR)))*fMRI_T
    # hrf=hrf[list(good_pts)]
    hrf=hrf[good_pts.astype(int)]
    # hrf = hrf([0:(p(7)/RT)]*fMRI_T + 1);
    hrf = hrf/np.sum(hrf);
    return hrf


if __name__ == "__main__":
    parser = argparse.ArgumentParser() 
    parser.add_argument('-fsl',dest='fsl',action='store_true',
            default=False,help='Whether it is an FSL-preprocessed dataset. Assumed true if fsldir is provided')
    parser.add_argument('--fsldir',dest='fsldir',default='',help='If FSL, path to FEAT Directory')
    parser.add_argument('--whichevs', dest='whichevs',type=int, nargs='*',help='If FSL, list of EVs to Compute Beta Series for. Number corresponds to real EV number in FEAT.')
    parser.add_argument('--numrealev',dest='numrealev',type=int,help='If FSL, total number of real EVs in FEAT model')
    parser.add_argument('-motpars',dest='motpars',action='store_true',
            default=False,help='If FSL, include motion parameters in model. The code assumes that motion parameters are the first 6 EVs (12 if including temporal derivative EVs) after the real EVs in the FEAT design matrix.')
    parser.add_argument('-tempderiv',dest='tempderiv',action='store_true',
            default=False,help='If FSL, include temporal derivates. The code assumes that temporal derivatives are immediately after each EV/motion parameter in the FEAT design matrix.')
    parser.add_argument('--image',dest='img',
            default='',help='If not using FSL, path to the preprocessed fMRI image')
    parser.add_argument('--mask', dest='mask',
            default='',help='If not using FSL, path to the brain mask for preprocessed image')
    parser.add_argument('--tr',dest='TR',type=float,
            default=None,help='If not using FSL, repetition time (TR) in seconds')
    parser.add_argument('--hpf',dest='HPF',type=float,
            default=None,help='If not using FSL, the highpass filter used in preprocessing, in seconds')
    parser.add_argument('--design',dest='design',
            default=None,help='If not using FSL, path to a csv file with conditions, onsets, durations (and amplitude if you want)')
    parser.add_argument('--confounds',dest='confounds',
            default=None,help='If not using FSL, path to a delimited file with confounds (e.g., motion parameters)')
    parser.add_argument('--outdir',dest='outdir',
            default='',help='If not using FSL, provide an output directory to save everything. Otherwise saved in a "betaseries" folder in the FEAT directory.')
    parser.add_argument('--timeres',dest='timeres',type=float,default=0.001, help='Optional- time resolution for convolution. (default= .001)')
    parser.add_argument('-LSA',dest='LSA',action='store_true',
            default=False,help='Include tag to compute LSA.')
    parser.add_argument('-LSS',dest='LSS',action='store_true',
            default=False,help='Include tag to compute LSS.')     

    args = parser.parse_args()
        
    
    if not args.fsl and args.fsldir !='':
        args.fsl = True
        
    pybeta=pybetaseries(fsfdir = args.fsldir,tempderiv = args.tempderiv, motpars = args.motpars, time_res = args.timeres, img = args.img, mask=args.mask, design= args.design, fsl=args.fsl,outdir = args.outdir,TR=args.TR,HPF = args.HPF)
    
    if args.LSS:    
        if args.fsl:
            for ev in args.whichevs:
                pybeta.LSS([ev],args.numrealev)
        else:
                pybeta.LSS(None,None)
    if args.LSA:
        pybeta.LSA(args.whichevs,args.numrealev)
        