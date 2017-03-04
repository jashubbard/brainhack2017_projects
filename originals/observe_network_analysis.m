%% load all the timeseries

clear
clc

alldata_raw = readtable('~/Dropbox/data_analysis/observe_fmri/extracted_data/network/timeseries_network_all.txt','Delimiter','\t');
alldata_raw = sortrows(alldata_raw,{'subject','vol'});


%also add timeseries of left and righ NAcc
otherrois = readtable('~/Dropbox/data_analysis/observe_fmri/extracted_data/timeseries_all.txt','Delimiter','\t');

alldata_raw = innerjoin(alldata_raw,otherrois(:,{'subject','vol','l_nacc','r_nacc'}),'keys',{'subject','vol'});

alldata_raw.Properties.VariableNames{end-1} = 'roi__10_14__6_nii_gz'; %left NAcc
alldata_raw.Properties.VariableNames{end} = 'roi_10_12__8_nii_gz'; %right NAcc

allvars = alldata_raw.Properties.VariableNames(3:end);
%% get all the filenames and extract the actual coordinates from the names

filenames = strrep(allvars,'_nii_gz','.nii.gz')';

cleanvars = strrep(allvars,'_nii_gz','');
cleanvars = strrep(cleanvars,'__','_-');
cleanvars = strrep(cleanvars,'roi_','');

coords = cellfun(@(x) strsplit(x,'_'),cleanvars,'Uniform',false)';
coords = str2double(vertcat(coords{:}));

%% create a table with all the roi info

orig_order = [1:length(allvars)]'; %the order of the original dataset (alphabetical)

roi_info = table(orig_order,allvars',filenames,cleanvars',coords(:,1),coords(:,2),coords(:,3),'VariableNames',{'orig_order','varname','filename','cleanname','x','y','z'});

%this has the network information, but with TAL coordinates
coleetal = readtable('~/Dropbox/data_analysis/observe_fmri/coords from cole et al.xlsx','ReadVariableNames',true);


%this has the ROI number, but with MNI coordinates
poweretal = readtable('~/Dropbox/data_analysis/observe_fmri/power_2011_coords.xlsx','ReadVariableNames',true);
poweretal.x = round(poweretal.x);
poweretal.y = round(poweretal.y);
poweretal.z = round(poweretal.z);

%fix this 1 coordinate to match up with our data
poweretal.y(poweretal.x==-2 & poweretal.y == -13 & poweretal.z== 12) = -12;


%get the ROI number from the Power et al. list
roi_info = outerjoin(roi_info,poweretal,'keys',{'x','y','z'},'MergeKeys',true);

%add the network column from the Cole et al. list
roi_info = outerjoin(roi_info,coleetal(:,{'ROI','Net'}),'keys','ROI','MergeKeys',true);

%make sure it's sorted in the original order, so we can reshuffle into the
%correct order
roi_info = sortrows(roi_info,'orig_order');

%need to add ROI and Network info for the NAcc rois
roi_info.Net(isnan(roi_info.ROI)) = 10; %subcortical network 
roi_info.ROI(isnan(roi_info.ROI)) = roi_info.orig_order(isnan(roi_info.ROI));


roi_info.Net_name = cell(height(roi_info),1);
roi_info.Net_name(roi_info.Net==-1) = {'unknown'};
roi_info.Net_name(roi_info.Net==1) = {'somato-motor'};
roi_info.Net_name(roi_info.Net==3) = {'cingulo-opercular'};
roi_info.Net_name(roi_info.Net==4) = {'auditory'};
roi_info.Net_name(roi_info.Net==5) = {'default-mode'};
roi_info.Net_name(roi_info.Net==7) = {'visual'};
roi_info.Net_name(roi_info.Net==8) = {'fronto-parietal'};
roi_info.Net_name(roi_info.Net==9) = {'salience'};
roi_info.Net_name(roi_info.Net==10) = {'subcortical'};
roi_info.Net_name(roi_info.Net==11) = {'ventral-attention'};
roi_info.Net_name(roi_info.Net==12) = {'dorsal-attention'};

%% reorder the columns into the order we want


%sort by roi # (from cole et al and power et al)
[roi_sorted,idx] = sortrows(roi_info,{'ROI'});

%shuffle the columns of alldata so they're in the same order as the ROI
%number
temp = alldata_raw(:,3:end);
temp2 = temp(:,idx);
alldata = horzcat(alldata_raw(:,1:2),temp2);



% 
% %% checking conversion from tal to mni (because same network used in different papers, but different spaces)
% poweretal = readtable('~/Dropbox/Observed_Altruism/braindata/power_2011_coords.xlsx','ReadVariableNames',true);
% 
% talpoints = icbm_other2tal(poweretal{:,{'x','y','z'}})
% 
% %%

%% read in design information

alldesign = readtable('~/Dropbox/data_analysis/observe_fmri/all_designmats.txt','Delimiter','\t');
alldesign = sortrows(alldesign,{'subject','vol'});

% %% get condition_specific timeseries (things get too correlated this way)
% 
% 
% sc1_charity = bsxfun(@times,alldata{:,3:end},alldesign{:,'sc1_charity'});
% sc1_win = bsxfun(@times,alldata{:,3:end},alldesign{:,'sc1_win'});
% sc1_lose = bsxfun(@times,alldata{:,3:end},alldesign{:,'sc1_lose'});
% sc1_baseline = bsxfun(@times,alldata{:,3:end},alldesign{:,'sc1_baseline'});
% 
% 
% sc2_charity = bsxfun(@times,alldata{:,3:end},alldesign{:,'sc2_charity'});
% sc2_win = bsxfun(@times,alldata{:,3:end},alldesign{:,'sc2_win'});
% sc2_lose = bsxfun(@times,alldata{:,3:end},alldesign{:,'sc2_lose'});
% sc2_baseline = bsxfun(@times,alldata{:,3:end},alldesign{:,'sc2_baseline'});





%% get correlation matricies for each condition

allsubs = unique(alldata.subject);

corr_charity1 = nan(length(allvars),length(allvars),length(allsubs));
corr_win1 = nan(size(corr_charity1));
corr_lose1 = nan(size(corr_charity1));
corr_baseline1 = nan(size(corr_charity1));

corr_charity2 = nan(size(corr_charity1));
corr_win2 = nan(size(corr_win1));
corr_lose2 = nan(size(corr_lose1));
corr_baseline2 = nan(size(corr_charity1));

for s=1:length(allsubs)
   
    rawdata = alldata{alldata.subject==allsubs(s),3:end};
    desdata = alldesign(alldata.subject==allsubs(s),:);
    
    %screen 1 conditions, get corrlation matricies for each
    sc1_charity = rawdata(desdata.sc1_charity>0,:);
    corr_charity1(:,:,s) = corr(sc1_charity);
    
    sc1_win = rawdata(desdata.sc1_win>0,:);
    corr_win1(:,:,s) = corr(sc1_win);
    
    sc1_lose = rawdata(desdata.sc1_lose>0,:);
    corr_lose1(:,:,s) = corr(sc1_lose);
    
    sc1_baseline = rawdata(desdata.sc1_baseline>0,:);
    corr_baseline1(:,:,s) = corr(sc1_baseline);
    
    %screen 2
    sc2_charity = rawdata(desdata.sc2_charity>0,:);
    corr_charity2(:,:,s) = corr(sc2_charity);
    
    sc2_win = rawdata(desdata.sc2_win>0,:);
    corr_win2(:,:,s) = corr(sc2_win);
    
    sc2_lose = rawdata(desdata.sc2_lose>0,:);
    corr_lose2(:,:,s) = corr(sc2_lose);
    
    sc2_baseline = rawdata(desdata.sc2_baseline>0,:);
    corr_baseline2(:,:,s) = corr(sc2_baseline);
    
    
    
      
end

%%

sub = randsample(size(corr_charity1,3),1);
close all

figure;

subplot(2,2,1)
imagesc(corr_charity1(:,:,sub));
colorbar;
title('charity')

subplot(2,2,2)
imagesc(corr_win1(:,:,sub));
colorbar;
title('win')


subplot(2,2,3)
imagesc(corr_lose1(:,:,sub));
colorbar;
title('lose')

subplot(2,2,4)
imagesc(corr_baseline1(:,:,sub));
colorbar;
title('baseline')


figure;

subplot(2,2,1)
imagesc(corr_charity2(:,:,sub));
colorbar;
title('charity2')

placeFigures

subplot(2,2,2)
imagesc(corr_win2(:,:,sub));
colorbar;
title('win2')


subplot(2,2,3)
imagesc(corr_lose2(:,:,sub));
colorbar;
title('lose2')

subplot(2,2,4)
imagesc(corr_baseline2(:,:,sub));
colorbar;
title('baseline2')


%% threshold and binarize

net_threshold = .10;

allbin = struct;
allbin.charity1 = nan(size(corr_charity1));
allbin.win1 = nan(size(allbin.charity1));
allbin.lose1 = nan(size(allbin.charity1));
allbin.baseline1 = nan(size(allbin.charity1));

allbin.charity2 = nan(size(allbin.charity1));
allbin.win2 = nan(size(allbin.charity1));
allbin.lose2 = nan(size(allbin.charity1));
allbin.baseline2 = nan(size(allbin.charity1));

for s=1:length(allsubs)
    
    %screen1
    thresh = threshold_proportional(corr_charity1(:,:,s),net_threshold);
    allbin.charity1(:,:,s) = weight_conversion(thresh, 'binarize');
    
    thresh = threshold_proportional(corr_win1(:,:,s),net_threshold);
    allbin.win1(:,:,s) = weight_conversion(thresh, 'binarize');
    
    thresh = threshold_proportional(corr_lose1(:,:,s),net_threshold);
    allbin.lose1(:,:,s) = weight_conversion(thresh, 'binarize');
    
    thresh = threshold_proportional(corr_baseline1(:,:,s),net_threshold);
    allbin.baseline1(:,:,s) = weight_conversion(thresh, 'binarize');
    
    %screen2
    thresh = threshold_proportional(corr_charity2(:,:,s),net_threshold);
    allbin.charity2(:,:,s) = weight_conversion(thresh, 'binarize');
    
    thresh = threshold_proportional(corr_win2(:,:,s),net_threshold);
    allbin.win2(:,:,s) = weight_conversion(thresh, 'binarize');
    
    thresh = threshold_proportional(corr_lose2(:,:,s),net_threshold);
    allbin.lose2(:,:,s) = weight_conversion(thresh, 'binarize');
    
    thresh = threshold_proportional(corr_baseline2(:,:,s),net_threshold);
    allbin.baseline2(:,:,s) = weight_conversion(thresh, 'binarize');
    
   
end

%%

allconds = {'charity1','win1','lose1','baseline1','charity2','win2','lose2','baseline2'};

numsubs = length(allsubs);
numnodes = size(corr_baseline1,2);


results = struct;
results.subs = allsubs;
results.roi_info = roi_sorted;

for c = 1:length(allconds)
    
    ds = allbin.(allconds{c});
    
    
    betweenness = nan(numsubs,numnodes);
    community_struct = nan(numsubs,numnodes);
    clustering = nan(numsubs,numnodes);
    modularity = nan(numsubs,1);
    participation = nan(numsubs,numnodes);
    Eglob = nan(numsubs,1);
    Eloc = nan(numsubs,numnodes);
    assort = nan(numsubs,1);
    GVC = nan(numsubs,numnodes);
    
    
    
    
    for i = 1:numsubs
        
        fprintf('...starting %s, subject %d\n',allconds{c},i);
        
        W_bin = ds(:,:,i);
        

        betweenness(i,:) = betweenness_bin(W_bin);
        %modularity
        [community_struct(i,:), modularity(i)] = modularity_louvain_und(W_bin);
        
        %participation coef
        participation(i,:) = participation_coef(W_bin,community_struct(i,:));
        
        %clustering coefficient
        clustering(i,:)=clustering_coef_bu(W_bin);
        
        %global and local efficiency
        Eglob(i) = efficiency_bin(W_bin);
        Eloc(i,:) = efficiency_bin(W_bin,1);
        
        %assortativity
        assort(i) = assortativity_bin(W_bin,0);
        
        %gvc value-- need original, unthresholded connectivity matrix
        conmat = eval(['corr_',allconds{c}]);
        conmat = conmat(:,:,i);
        GVC(i,:) = gvc(conmat);

        
        
    end
    
%     GVC(:,37) = []; %bad node
    
    %store all the results into a structure
    results.(allconds{c}).betweenness = betweenness;
    results.(allconds{c}).community_struct = community_struct;
    results.(allconds{c}).modularity = modularity;
    results.(allconds{c}).participation = participation;
    results.(allconds{c}).clustering = clustering;
    results.(allconds{c}).eglob = Eglob;
    results.(allconds{c}).eloc = Eloc;
    results.(allconds{c}).assort = assort;
    results.(allconds{c}).gvc = GVC;
    results.(allconds{c}).wbin = ds;
   
    fprintf('finished %s\n',allconds{c});
    
end

%%

results.coords = roi_sorted{:,{'x','y','z'}};

cd ~/Dropbox/data_analysis/observe_fmri/
save('network_results.mat','-struct','results')


%%

clear
cd ~/Dropbox/data_analysis/observe_fmri/
load('network_results.mat')

results = v2struct;

behdata = readtable('~/Dropbox/data_analysis/observe_fmri/all_observe_data.csv','Delimiter',',');

beh = behdata(:,{'subject','age','rate_diff','knorms','Ugiving','Ogiving'});

beh.glob_charity1 = results.charity1.eglob;
beh.glob_win1 = results.win1.eglob;

beh.glob_charity2 = results.charity2.eglob;
beh.glob_win2 = results.win2.eglob;


beh.mod_charity1 = results.charity1.modularity;
beh.mod_win1 = results.win1.modularity;

beh.mod_charity2 = results.charity2.modularity;
beh.mod_win2 = results.win2.modularity;
%%

roi_sorted = results.roi_info;

for i = 1:2
    
   [~,subcort_res] = system(sprintf('source ~/.bash_profile; atlasquery -a "Harvard-Oxford Subcortical Structural Atlas" -c %d,%d,%d',roi_sorted.x(i),roi_sorted.y(i),roi_sorted.z(i)))
   [~,cort_res] = system(sprintf('atlasquery -a "Harvard-Oxford Cortical Structural Atlas" -c %d,%d,%d',roi_sorted.x(i),roi_sorted.y(i),roi_sorted.z(i)))
    
    
end




%% correlations

behvar = 'Ugiving'
condition = 'charity1'
measure = 'participation'


cors = table;
cors.r = nan(numnodes,1);
cors.p = nan(numnodes,1);


for i = 1:numnodes
   
    [r,p] = corrcoef(beh.(behvar),results.(condition).(measure)(:,i));
    
    cors.r(i) = r(1,2);
    cors.p(i) = p(1,2);
    
    
    
end


find(cors.p<=.06)
results.roi_info.Net_name(cors.p<=.06)


