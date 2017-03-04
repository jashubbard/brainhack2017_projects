# -*- coding: utf-8 -*-
"""
Created on Wed Oct 15 11:28:24 2014

@author: jason
"""
#%%

import networkx as nx
import numpy as np
from mayavi import mlab
import scipy.io.matlab as sio
import os


#load data
os.chdir(os.path.expanduser('~/Dropbox/data_analysis/observe_fmri'))
temp=sio.loadmat('network_results.mat')


condition = 'lose1'
allnets=['subcortical','cingulo-opercular','fronto-parietal','dorsal-attention','ventral-attention','auditory','salience']
viewnetwork = 'all'
#viewnetwork = 'all'
group = True

subnum=10


measure='betweenness'
#measure='eloc'
#measure='participation'
#measure='clustering'
#measure='community_struct'




bigmatrix = temp[condition]['wbin'][0][0]


#%%



if group:
    H=nx.to_networkx_graph(bigmatrix[:,:,0])
else:
    H=nx.to_networkx_graph(bigmatrix[:,:,subnum])
    


# reorder nodes from 0,len(G)-1
G=nx.convert_node_labels_to_integers(H,first_label=0)
# 3d spring layout
# pos=nx.spring_layout(G,dim=3)
# numpy array of x,y,z positions in sorted node order
# xyz=np.array([pos[v] for v in sorted(G)])
# scalar colors
# scalars=np.array(G.nodes())+5


#calculating stuff
#clustering
#scalars=np.power(np.array(nx.clustering(G).values()),2)

#degree
# scalars=np.array(nx.degree(G).values())

#betweenness centrality
scalars=np.array(nx.betweenness_centrality(G).values())*100







scale_factors={'participation': 15.0,
               'eloc': 10,
               'betweenness': .01,
               'community_struct': .5,
               'clustering': 15}


if group:
#average of group
    vals = np.mean(temp[condition][measure][0][0],axis=0)
    
    scalars = np.zeros(len(scalars))
    
    if viewnetwork !='all':
        idx = np.array([x in viewnetwork for x in temp['net_name']])
        scalars[idx] = vals[idx]
    else:
        scalars = vals
        
else:
    #data from a single subject
    scalars=np.array(temp[condition][measure][0][0][subnum,:],dtype='float');






#xyz=np.loadtxt("AAL_coords.txt", delimiter="\t")
xyz = temp['coords']
# xyz=xyz[40:len(xyz),:]


mlab.close(all=True)
mlab.figure(bgcolor=(0, 0, 0))
# mlab.clf()

pts = mlab.points3d(xyz[:,0], xyz[:,1], xyz[:,2],
                    scalars,
                    scale_factor= scale_factors[measure],
                    scale_mode='scalar',
                    colormap='jet',
                    resolution=20)

pts.module_manager.scalar_lut_manager.show_scalar_bar = True


pts.mlab_source.dataset.lines = np.array(G.edges())
tube = mlab.pipeline.tube(pts, tube_radius=0.3)
mlab.pipeline.surface(tube, color=(0.5, 0.5, 0.5),opacity=0.1)


#mlab.savefig('mayavi2_spring.png')
mlab.show() ##interactive window
