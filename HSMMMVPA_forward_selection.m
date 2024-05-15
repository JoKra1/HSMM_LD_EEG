% Estimate HSMMs for Lexical Decision
clear
clc
close all

subs = 1:29;

% Optional exclusion.
exclude = [6,15,19,20,26];
subs = subs(~ismember(subs,exclude));
nsubs = numel(subs);

% Paths assume that repository is working directory/current folder.

% Path to analysis scripts, also used to write results
analysis_path = './';

% Path to data
saveTo = './data/';

cd(analysis_path);

%Add HSMM functions to path (Download from: https://osf.io/z49me/files/)
if ~contains(path,'HsMM_functions')
  addpath('HsMM_functions');
end

% Specify data to load
comp_add = '_excl_6_15_bs_100_';

% EEGlab is required to plot the recovered stage topologies
% can be downloaded at: https://eeglab.org/others/How_to_download_EEGLAB.html
if ~contains(path,'eeglab') % Might require re-naming
  addpath('eeglab');
end

eeglab;

%% PCA CHECK
close all
% load data
load([saveTo strcat('data4HMM_baseline_LD',comp_add,'.mat')]);

plot_path = strcat('results/plots/');
mkdir([analysis_path plot_path]);

%Check variance explained in PCAs
tmp = cumsum(latent10/sum(latent10)*100);
tmp(1:20)'
plot(cumsum(latent10/sum(latent10)*100))
saveas(gcf,strcat(analysis_path,plot_path,'pca_var.png'))
close
% with 10 there is >= 90% variance explained

%% Find initial parameters for each fold
mkdir([analysis_path strcat('results/models/')]);
mkdir([analysis_path strcat('results/init_params/')]);
res_path = strcat('results/init_params/LD_baseline_all_sols_singleModel/');

% res_path holds initial parameters for every fold.
mkdir([analysis_path res_path]);

pcs = 10; %number of PCs to use.
c=0; %no correlation constraints

parpool(5)

%select data
score = normedscore10(:,1:pcs);

%rebuild x/y
trial_lens = y-x+1; %calc lengths

max_dur = max(trial_lens);   
mean_dur = mean(trial_lens);
min_dur = min(trial_lens);
max_bumps = floor(min_dur / 5);
disp(min_dur)
disp(max_bumps)

% Collect starting parameters for every fold
optimal_bumps = {nsubs,max_bumps};
optimal_gammas = {nsubs,max_bumps};

%initialize gamma params (equally spaced bumps on every trial)
for Nbump = 1:max_bumps
    params8{Nbump} = repmat([2 ceil(max_dur)/(Nbump+1)/2], Nbump+1,1); %set initial gamma params to max RT / nr of stages
end

%fit models from max bumps to 1 bump, to find best initial params for
%LOOCV for every fold separately

parfor subj = 1:nsubs
    %select data rest - i.e., without subj
    score_rest = score(subjects_varS ~= subj,:);

    %x/y rest
    trial_lens_rest = trial_lens(subjects_var ~= subj); 
    x_rest = [0; cumsum(trial_lens_rest)] + 1;
    x_rest(end) = [];
    y_rest = x_rest + trial_lens_rest - 1;

    % Back-fit model
    [lkhs, mags, params, lkhsb, magsb,...
    paramsb, lkhsI, magsI, paramsI, lkhsU, magsU,...
    paramsU, Is] = fitBumpsIU_m(score_rest, x_rest, y_rest, params8, max_bumps);

    %for each solution...
    for Nbumps = 1:max_bumps
    
        %get best params
        best_bumps = magsb{Nbumps};
        best_gammas = paramsb{Nbumps};
    
        optimal_bumps{subj,Nbumps} = best_bumps;
        optimal_gammas{subj,Nbumps} = best_gammas;
    end
end

% Save init parameters for all folds
save(strcat(analysis_path, res_path, 'optimal_sols_baseline_singleModel',comp_add,'.mat'), 'optimal_bumps', 'optimal_gammas');

delete(gcp('nocreate'))


%% LOOCV - based on optimal solutions above
%To obtain the optimal number of bumps we apply LOOCV (make HSMM model for all subjects...
% but one and test it on this one participant. We apply that in a for loop to each participant).

parpool(5)

% shared model solution
res_path = strcat('results/init_params/LD_baseline_all_sols_singleModel/');
load(strcat(analysis_path, res_path, 'optimal_sols_baseline_singleModel',comp_add,'.mat'));

%select data
score = normedscore10(:,1:pcs);

%rebuild x/y
trial_lens = y-x+1; %calc lengths

max_dur = max(trial_lens);
min_dur = min(trial_lens);
max_bumps = floor(min_dur / 5);

for Nbumps = 1:max_bumps

    disp([' - Nbumps ', num2str(Nbumps)])

    likelihoods = repmat(-Inf,nsubs,1);
    bumps = zeros(pcs,Nbumps,nsubs);
    gammas = zeros(Nbumps+1, 2, nsubs);

    parfor subj = 1:nsubs

        %select data one (the data from the left-out subject)
        score_one = score(subjects_varS == subj,:);

        %x/y one
        trial_lens_one = trial_lens(subjects_var == subj);
        x_one = [0; cumsum(trial_lens_one)] + 1;
        x_one(end) = [];
        y_one = x_one + trial_lens_one - 1;

        % Take the fold-specific back-fitting estimates, no need
        % to fit the models again here.
        bumps(:,:,subj) = optimal_bumps{subj,Nbumps};
        gammas(:,:,subj) = optimal_gammas{subj,Nbumps};

        % evaluate the model on the left out subject
        [likelihoods(subj), ~, ~, ~] = hsmmBoundCorrsVary(score_one,bumps(:,:,subj), gammas(:,:,subj),0,x_one,y_one,max_dur,c);
    end

    % calc final solution based on average
    % subject
    [~,bumps_final,gammas_final, events_final]= hsmmBoundCorrsVary(score,nanmean(bumps,3),nanmean(gammas,3),1,x,y,max_dur,c);

    % save to
    save(strcat(analysis_path, 'results/models/loocv_baseline_singleModel_bumps', num2str(Nbumps),comp_add,'.mat'),...
   'likelihoods', 'bumps', 'gammas', 'bumps_final', 'gammas_final','events_final');

end

delete(gcp('nocreate'))

%Write trial-level results
pcs=10;
ncond=3;

%select data
score = normedscore10(:,1:pcs);

%rebuild x/y
trial_lens = y-x+1; %calc lengths

max_dur = max(trial_lens);  % duration of each trial
min_dur = min(trial_lens);
max_bumps = floor(min_dur / 5);


loocv = {};
likelis = zeros(nsubs,max_bumps);
for Nbumps = 1:max_bumps
    loocv{Nbumps} = load(strcat(analysis_path, 'results/models/loocv_baseline_singleModel_bumps', num2str(Nbumps),comp_add,'.mat'));
    likelis(:,Nbumps) = loocv{Nbumps}.likelihoods;
end

%improvement
improv = zeros(Nbumps,1);
for Nbumps = 2:max_bumps
    improv(Nbumps) = sum(likelis(:,Nbumps) > likelis(:,Nbumps-1));
end
improv = [(1:Nbumps)' improv];

disp(improv)

for Nbumps = 1:max_bumps
    
    tmp = loocv{Nbumps};
       
    %calc durations based on event probs
    
    %first, calc eventProbs
    best_events = tmp.events_final; %samples x trials x params
    
    nTrials = size(best_events,2);
    maxDur = size(best_events,1);
    
    %first get onsets + response
    onsets = zeros(nTrials, Nbumps); %trials x bumps
    for tr = 1:nTrials
        onsets(tr,:) = (1:maxDur) * squeeze(best_events(:,tr,:));
    end
    onsets(:,end+1) = y-x; % add end trials
    onsets = onsets .* 10; % got to ms
    
    %calc stage durations
    durations = onsets;
    for b = 2:(Nbumps+1)
        durations(:,b) = onsets(:,b) - onsets(:,b-1);
    end
    
    % write results to file for visualization in R
    trial_dur_exp = table(durations,onsets,subjects_var,trials_var,rts_var,blocks_var);
    writetable(trial_dur_exp,strcat(analysis_path,'results/trial_hsmm_dat_single_model_',num2str(Nbumps),comp_add,'.csv'));
end

%% Likelihood comparison I

% to save plots
res_path = strcat('results/plots/');

pcs = 10; %number of PCs to use.
c=0; %no correlation constraints
max_bumps = 6;

% Collect likilhood scores
loocv = {};
likelis_overall = zeros(nsubs,max_bumps);
for Nbumps = 1:max_bumps
    loocv{Nbumps} = load(strcat(analysis_path, 'results/models/loocv_baseline_singleModel_bumps', num2str(Nbumps),comp_add,'.mat'));
    likelis_overall(:,Nbumps) = loocv{Nbumps}.likelihoods;
end

%improvement
improv_overall = zeros(Nbumps,1);
for Nbumps = 2:max_bumps
    Nbumps
    improv_overall(Nbumps) = sum(likelis_overall(:,Nbumps) > likelis_overall(:,Nbumps-1));
    [p, h] = signtest(likelis_overall(:,Nbumps),likelis_overall(:,Nbumps-1),'Tail','right')
end
nbumps = (1:Nbumps)';
heldoutllk = mean(likelis_overall,1)';
llk_table = table(nbumps,improv_overall,heldoutllk)
writetable(llk_table,strcat(analysis_path,'results/llk_table_',comp_add,'.csv'));
improv_overall = [nbumps improv_overall];

%% likelihood plot I

% Over-all
like_to_plot = 6;

hold on

xlabel('Number of Bumps','fontsize',18)
xticks(1:(Nbumps));
xlim([.5,like_to_plot+.5]);
ylabel('Log-likelihood','fontsize',18)

pl = plot(mean(likelis_overall(:,1:like_to_plot)),'.-');

for b = 2:max_bumps
    text(pl.XData(b),pl.YData(b)+.05*range(pl.YData),strcat(num2str(improv_overall(b,2)),'/',num2str(nsubs)),'FontSize',18,'HorizontalAlignment','center','Color', pl.Color)
end
pl.MarkerSize = 32;

overall_col = pl.Color;

ylim([-3000,-1300]);

set(gcf,'Color',[1 1 1]);
set(gcf,'Units','Inches');
pos = get(gcf,'Position');
set(gcf, 'PaperPositionMode', 'auto','PaperUnits','inches','PaperSize',[pos(3), pos(4)]);
print('-dpdf', strcat(analysis_path, res_path, 'likelihoods_comparison_overall',comp_add,'.pdf'));
close;

%% First forward pass
% Gamma parameters differ between conditions for every stage.
map_I = [1 1 1 1 1 1;...
         2 1 1 1 1 1;...
         3 1 1 1 1 1];

map_II = [1 1 1 1 1 1;...
          1 2 1 1 1 1;...
          1 3 1 1 1 1];

map_III = [1 1 1 1 1 1;...
           1 1 2 1 1 1;...
           1 1 3 1 1 1];

map_IV = [1 1 1 1 1 1;...
          1 1 1 2 1 1;...
          1 1 1 3 1 1];

map_V = [1 1 1 1 1 1;...
         1 1 1 1 2 1;...
         1 1 1 1 3 1];

map_VI = [1 1 1 1 1 1;...
          1 1 1 1 1 2;...
          1 1 1 1 1 3];

maps = {map_I,map_II,map_III,map_IV,map_V,map_VI};
map_indices = ["I","II","III","IV","V","VI"];

% re-load data
load([saveTo strcat('data4HMM_baseline_LD',comp_add,'.mat')]);

%load optimal solutions from fitbumps overall model
res_path = strcat('results/init_params/LD_baseline_all_sols_singleModel/');
load(strcat(analysis_path, res_path, 'optimal_sols_baseline_singleModel',comp_add,'.mat'));

%select data
pcs=10;
c=0;
score = normedscore10(:,1:pcs);

%rebuild x/y
trial_lens = y-x+1; %calc lengths //dur

max_dur = max(trial_lens);
min_dur = min(trial_lens);
max_bumps = floor(min_dur / 5);

% Number of word types
ncond = 3;

% Get likelihood from current best model:
loocv_best = load(strcat(analysis_path, 'results/models/loocv_baseline_singleModel_bumps5',comp_add,'.mat'));
likelis_best = loocv_best.likelihoods;

for map_i = 1:length(maps)
    map_index = map_indices(map_i)
    cur_map = maps{map_i};

    %specify nr of bumps
    Nbumps = sum(cur_map(1,:))-1;
    
    disp(['Map - Nbumps ', num2str(Nbumps)])
    
    likelihoods = repmat(-Inf,nsubs,1); %BestLikelihood
    bumps = zeros(pcs,size(cur_map,2)-1,nsubs); %BestBumpMag
    gammas = zeros(size(cur_map,2), 2, size(cur_map,1), nsubs); %BestGamma
    
    parpool(5)
    
    parfor subj = 1:nsubs
        
        % Select starting parameters for this fold.
        new_bumps = optimal_bumps{subj,Nbumps};
    
        new_gammas = optimal_gammas{subj,Nbumps};
    
        %select data one
        score_one = score(subjects_varS == subj,:);
    
        %x/y one
        trial_lens_one = trial_lens(subjects_var == subj);
        x_one = [0; cumsum(trial_lens_one)] + 1;
        x_one(end) = [];
        y_one = x_one + trial_lens_one - 1;
        conds_one = conditions(subjects_var == subj);
        
    
        %select data rest
        score_rest = score(subjects_varS ~= subj,:);
    
        %x/y rest
        trial_lens_rest = trial_lens(subjects_var ~= subj); 
        x_rest = [0; cumsum(trial_lens_rest)] + 1;
        x_rest(end) = [];
        y_rest = x_rest + trial_lens_rest - 1;
        conds_rest = conditions(subjects_var ~= subj);
    
        %calculate model based on other subjects
        [~, bumps(:,:,subj), gammas(:,:,:,subj),~]=hsmmBoundCorrsCondsVary(score_rest, conds_rest, ...
                                                                           new_bumps,new_gammas, ...
                                                                           1,x_rest,y_rest,max_dur,cur_map,c);
            
    
        % fit the model on the left out subject
        [likelihoods(subj), ~, ~, ~] = hsmmBoundCorrsCondsVary(score_one,conds_one, bumps(:,:,subj), gammas(:,:,:,subj),0,x_one,y_one,max_dur, cur_map, c);
    end
    
    % calc final solution based on average
    % subject
    [~,bumps_final,gammas_final, events_final]= hsmmBoundCorrsCondsVary(score,conditions,nanmean(bumps,3),nanmean(gammas,4),1,x,y,max_dur,cur_map,c);
    
    % save to
    save(strcat(analysis_path, 'results/models/loocv_baseline_map_',map_index,'_bumps', num2str(Nbumps),comp_add,'.mat'),...
    'likelihoods', 'bumps', 'gammas', 'bumps_final', 'gammas_final','events_final');
    
    % Check improvement according to current candidate
    improv = sum(likelihoods > likelis_best)
    [p, h] = signtest(likelihoods,likelis_best,'Tail','right')
    
    delete(gcp('nocreate'))


end

%% Second forward pass
% Current best model has different Gamma parameter for stage five.
% For all stages except five we again check whether different Gamma
% parameters are beneficial. For the fifth stage we check whether
% word type specific bumps preceding stage six are beneficial.

map_I = [1 1 1 1 1 1;...
         2 1 1 1 2 1;...
         3 1 1 1 3 1];

map_II = [1 1 1 1 1 1;...
          1 2 1 1 2 1;...
          1 3 1 1 3 1];

map_III = [1 1 1 1 1 1;...
           1 1 2 1 2 1;...
           1 1 3 1 3 1];

map_IV = [1 1 1 1 1 1;...
          1 1 1 2 2 1;...
          1 1 1 3 3 1];

map_V = [1 1 1 1 1 0 0 1;...
         1 1 1 1 0 1 0 1;...
         1 1 1 1 0 0 1 1];

map_VI = [1 1 1 1 1 1;...
          1 1 1 1 2 2;...
          1 1 1 1 3 3];

maps = {map_I,map_II,map_III,map_IV,map_V,map_VI};
map_indices = ["Ia","IIa","IIIa","IVa","Va","VIa"];

% re-load data
load([saveTo strcat('data4HMM_baseline_LD',comp_add,'.mat')]);

%load optimal solutions from fitbumps overall model
res_path = strcat('results/init_params/LD_baseline_all_sols_singleModel/');
load(strcat(analysis_path, res_path, 'optimal_sols_baseline_singleModel',comp_add,'.mat'));

%select data
pcs=10;
c=0;
score = normedscore10(:,1:pcs);

%rebuild x/y
trial_lens = y-x+1; %calc lengths //dur

max_dur = max(trial_lens);
min_dur = min(trial_lens);
max_bumps = floor(min_dur / 5);

% Number of word types
ncond = 3;

% Get likelihood from current best model:
loocv_best = load(strcat(analysis_path, 'results/models/loocv_baseline_map_V_bumps5',comp_add,'.mat'));
likelis_best = loocv_best.likelihoods;

for map_i = 1:length(maps)
    map_index = map_indices(map_i)
    cur_map = maps{map_i};

    %specify nr of bumps
    Nbumps = sum(cur_map(1,:))-1;
    
    disp(['Map - Nbumps ', num2str(Nbumps)])
    
    likelihoods = repmat(-Inf,nsubs,1); %BestLikelihood
    bumps = zeros(pcs,size(cur_map,2)-1,nsubs); %BestBumpMag
    gammas = zeros(size(cur_map,2), 2, size(cur_map,1), nsubs); %BestGamma
    
    parpool(5)
    
    parfor subj = 1:nsubs
        
        % Select starting parameters for this fold
        new_bumps = optimal_bumps{subj,Nbumps};
        new_gammas = optimal_gammas{subj,Nbumps};
    
        colsums = sum(cur_map,1);
        if ismember(1,colsums)
            have_bump = [];
            ci = 1;
            col = 1;
            % Identify which stages have a different bump per word type
            while ci < length(colsums)
                if colsums(ci) == 1
                    have_bump(end+1) = col;
                    ci = ci + (ncond -1);
                end
                col = col + 1;
                ci = ci + 1;
            end
            
            % Now add for each word type and identified stage extra bump and
            % gamma parameters.
            if ~ismember(1,have_bump)
                tmp_bumps = new_bumps(:,1);
                tmp_gammas = new_gammas(1,:);
            else
                tmp_bumps = repmat(new_bumps(:,1),1,ncond);
                tmp_gammas = repmat(new_gammas(1,:),ncond,1);
            end
            
            for stage_i = 2:(Nbumps+1)
                
                if ~ismember(stage_i,have_bump)
                    if stage_i <= Nbumps % One bump less than number of stages!
                        tmp_bumps = horzcat(tmp_bumps,new_bumps(:,stage_i));
                    end
                    tmp_gammas = vertcat(tmp_gammas,new_gammas(stage_i,:));
                else
                    if stage_i <= Nbumps
                        tmp_bumps = horzcat(tmp_bumps,repmat(new_bumps(:,stage_i),1,ncond));
                    end
                    tmp_gammas = vertcat(tmp_gammas,repmat(new_gammas(stage_i,:),ncond,1));
                end
            end
            new_bumps = tmp_bumps
            new_gammas = tmp_gammas
        end
    
        %select data one
        score_one = score(subjects_varS == subj,:);
    
        %x/y one
        trial_lens_one = trial_lens(subjects_var == subj);
        x_one = [0; cumsum(trial_lens_one)] + 1;
        x_one(end) = [];
        y_one = x_one + trial_lens_one - 1;
        conds_one = conditions(subjects_var == subj);
        
    
        %select data rest
        score_rest = score(subjects_varS ~= subj,:);
    
        %x/y rest
        trial_lens_rest = trial_lens(subjects_var ~= subj); 
        x_rest = [0; cumsum(trial_lens_rest)] + 1;
        x_rest(end) = [];
        y_rest = x_rest + trial_lens_rest - 1;
        conds_rest = conditions(subjects_var ~= subj);
    
        %calculate model based on other subjects
        [~, bumps(:,:,subj), gammas(:,:,:,subj),~]=hsmmBoundCorrsCondsVary(score_rest, conds_rest, ...
                new_bumps,new_gammas, ...
                1,x_rest,y_rest,max_dur,cur_map,c);
            
    
        % fit the model on the left out subject
        [likelihoods(subj), ~, ~, ~] = hsmmBoundCorrsCondsVary(score_one,conds_one, bumps(:,:,subj), gammas(:,:,:,subj),0,x_one,y_one,max_dur, cur_map, c);
    end
    
    % calc final solution based on average
    % subject
    [~,bumps_final,gammas_final, events_final]= hsmmBoundCorrsCondsVary(score,conditions,nanmean(bumps,3),nanmean(gammas,4),1,x,y,max_dur,cur_map,c);
    
    % save to
    save(strcat(analysis_path, 'results/models/loocv_baseline_map_',map_index,'_bumps', num2str(Nbumps),comp_add,'.mat'),...
    'likelihoods', 'bumps', 'gammas', 'bumps_final', 'gammas_final','events_final');
    
    % Check improvement according to current candidate
    improv = sum(likelihoods > likelis_best)
    [p, h] = signtest(likelihoods,likelis_best,'Tail','right')
    
    delete(gcp('nocreate'))


end

%% Third forward pass 
% Current best model has different Gamma parameters for stage four and
% five. For stages 1,2,3, and 6 we again check whether different Gamma
% parameters are beneficial. For stages four and five we check whether
% varying bumps preceding stages five and six are beneficial.

map_I = [1 1 1 1 1 1;...
         2 1 1 2 2 1;...
         3 1 1 3 3 1];
    
map_II = [1 1 1 1 1 1;...
          1 2 1 2 2 1;...
          1 3 1 3 3 1];

map_III = [1 1 1 1 1 1;...
           1 1 2 2 2 1;...
           1 1 3 3 3 1];

map_IV = [1 1 1 1 0 0 1 1;...
          1 1 1 0 1 0 2 1;...
          1 1 1 0 0 1 3 1];

map_V = [1 1 1 1 1 0 0 1;...
         1 1 1 2 0 1 0 1;...
         1 1 1 3 0 0 1 1];

map_VI = [1 1 1 1 1 1;...
          1 1 1 2 2 2;...
          1 1 1 3 3 3];

maps = {map_I,map_II,map_III,map_IV,map_V,map_VI};
map_indices = ["Ib","IIb","IIIb","IVb","Vb","VIb"];

% re-load data
load([saveTo strcat('data4HMM_baseline_LD',comp_add,'.mat')]);

%load optimal solutions from fitbumps overall model
res_path = strcat('results/init_params/LD_baseline_all_sols_singleModel/');
load(strcat(analysis_path, res_path, 'optimal_sols_baseline_singleModel',comp_add,'.mat'));

%select data
pcs=10;
c=0;
score = normedscore10(:,1:pcs);

%rebuild x/y
trial_lens = y-x+1; %calc lengths //dur

max_dur = max(trial_lens);
min_dur = min(trial_lens);
max_bumps = floor(min_dur / 5);

% Number of word types
ncond = 3;

% Get likelihood from current best model:
loocv_best = load(strcat(analysis_path, 'results/models/loocv_baseline_map_IVa_bumps5',comp_add,'.mat'));
likelis_best = loocv_best.likelihoods;

for map_i = 1:length(maps)
    map_index = map_indices(map_i)
    cur_map = maps{map_i};

    %specify nr of bumps
    Nbumps = sum(cur_map(1,:))-1;
    
    disp(['Map - Nbumps ', num2str(Nbumps)])
    
    likelihoods = repmat(-Inf,nsubs,1); %BestLikelihood
    bumps = zeros(pcs,size(cur_map,2)-1,nsubs); %BestBumpMag
    gammas = zeros(size(cur_map,2), 2, size(cur_map,1), nsubs); %BestGamma

    parpool(5)
    
    parfor subj = 1:nsubs

        % Select starting parameters for this fold
        new_bumps = optimal_bumps{subj,Nbumps};
        new_gammas = optimal_gammas{subj,Nbumps};
    
        colsums = sum(cur_map,1);
        if ismember(1,colsums)
            have_bump = [];
            ci = 1;
            col = 1;
            % Identify which stages have a different bump per word type
            while ci < length(colsums)
                if colsums(ci) == 1
                    have_bump(end+1) = col;
                    ci = ci + (ncond -1);
                end
                col = col + 1;
                ci = ci + 1;
            end
            
            % Now add for each word type and identified stage extra bump and
            % gamma parameters.
            if ~ismember(1,have_bump)
                tmp_bumps = new_bumps(:,1);
                tmp_gammas = new_gammas(1,:);
            else
                tmp_bumps = repmat(new_bumps(:,1),1,ncond);
                tmp_gammas = repmat(new_gammas(1,:),ncond,1);
            end
            
            for stage_i = 2:(Nbumps+1)
                
                if ~ismember(stage_i,have_bump)
                    if stage_i <= Nbumps % One bump less than number of stages!
                        tmp_bumps = horzcat(tmp_bumps,new_bumps(:,stage_i));
                    end
                    tmp_gammas = vertcat(tmp_gammas,new_gammas(stage_i,:));
                else
                    if stage_i <= Nbumps
                        tmp_bumps = horzcat(tmp_bumps,repmat(new_bumps(:,stage_i),1,ncond));
                    end
                    tmp_gammas = vertcat(tmp_gammas,repmat(new_gammas(stage_i,:),ncond,1));
                end
            end
            new_bumps = tmp_bumps
            new_gammas = tmp_gammas
        end
    
        %select data one
        score_one = score(subjects_varS == subj,:);
    
        %x/y one
        trial_lens_one = trial_lens(subjects_var == subj);
        x_one = [0; cumsum(trial_lens_one)] + 1;
        x_one(end) = [];
        y_one = x_one + trial_lens_one - 1;
        conds_one = conditions(subjects_var == subj);
        
    
        %select data rest
        score_rest = score(subjects_varS ~= subj,:);
    
        %x/y rest
        trial_lens_rest = trial_lens(subjects_var ~= subj); 
        x_rest = [0; cumsum(trial_lens_rest)] + 1;
        x_rest(end) = [];
        y_rest = x_rest + trial_lens_rest - 1;
        conds_rest = conditions(subjects_var ~= subj);
    
        %calculate model based on other subjects
        [~, bumps(:,:,subj), gammas(:,:,:,subj),~]=hsmmBoundCorrsCondsVary(score_rest, conds_rest, ...
                new_bumps,new_gammas, ...
                1,x_rest,y_rest,max_dur,cur_map,c);
            
    
        % fit the model on the left out subject
        [likelihoods(subj), ~, ~, ~] = hsmmBoundCorrsCondsVary(score_one,conds_one, bumps(:,:,subj), gammas(:,:,:,subj),0,x_one,y_one,max_dur, cur_map, c);
    end
    
    % calc final solution based on average
    % subject
    [~,bumps_final,gammas_final, events_final]= hsmmBoundCorrsCondsVary(score,conditions,nanmean(bumps,3),nanmean(gammas,4),1,x,y,max_dur,cur_map,c);
    
    % save to
    save(strcat(analysis_path, 'results/models/loocv_baseline_map_',map_index,'_bumps', num2str(Nbumps),comp_add,'.mat'),...
    'likelihoods', 'bumps', 'gammas', 'bumps_final', 'gammas_final','events_final');
    
    % Check improvement according to current candidate
    improv = sum(likelihoods > likelis_best)
    [p, h] = signtest(likelihoods,likelis_best,'Tail','right')
    
    delete(gcp('nocreate'))


end

%% Fourth forward pass
% Strong evidence in favor of bump preceding stage five to vary for different
% word-types.
% Now we check again for additional Gammas and additional bumps.

map_I = [1 1 1 1 0 0 1 1;...
         2 1 1 0 1 0 2 1;...
         3 1 1 0 0 1 3 1];

map_II = [1 1 1 1 0 0 1 1;...
          1 2 1 0 1 0 2 1;...
          1 3 1 0 0 1 3 1];

map_III = [1 1 1 1 0 0 1 1;...
           1 1 2 0 1 0 2 1;...
           1 1 3 0 0 1 3 1];

map_IV = [1 1 1 1 0 0 1 0 0 1;...
          1 1 1 0 1 0 0 1 0 1;...
          1 1 1 0 0 1 0 0 1 1];

map_V = [1 1 1 1 0 0 1 1;...
         1 1 1 0 1 0 2 2;...
         1 1 1 0 0 1 3 3];

% Additionally, we need to check whether an additional bump should be
% inserted before/after the fifth stage for any of the word-types. This is
% because we found the bump preceding stage five to differ between
% word-types, which implies that the event in the EEG marking the onset of
% stage five looks very different between word types. We now need to check
% whether "looks very different" means that the stage should be split/an
% extra stage should be inserted for any individual word type.

map_VI = [1 1 1 1 0 0 1 1 1;... % Extra stage for words
          1 1 1 0 1 0 0 2 1;...
          1 1 1 0 0 1 0 3 1];

map_VII = [1 1 1 1 0 0 0 1 1;... % Extra stage for pseudo-words
           1 1 1 0 1 0 1 2 1;...
           1 1 1 0 0 1 0 3 1];

map_VIII = [1 1 1 1 0 0 0 1 1;... % Extra stage for random strings
            1 1 1 0 1 0 0 2 1;...
            1 1 1 0 0 1 1 3 1];

% First test the simple maps...

maps = {map_I,map_II,map_III,map_IV,map_V};
map_indices = ["Ic","IIc","IIIc","IVc","Vc"];

% re-load data
load([saveTo strcat('data4HMM_baseline_LD',comp_add,'.mat')]);

%load optimal solutions from fitbumps overall model
res_path = strcat('results/init_params/LD_baseline_all_sols_singleModel/');
load(strcat(analysis_path, res_path, 'optimal_sols_baseline_singleModel',comp_add,'.mat'));

%select data
pcs=10;
c=0;
score = normedscore10(:,1:pcs);

%rebuild x/y
trial_lens = y-x+1; %calc lengths //dur

max_dur = max(trial_lens);
min_dur = min(trial_lens);
max_bumps = floor(min_dur / 5);

% Number of word types
ncond = 3;

% Get likelihood from current best model:
loocv_best = load(strcat(analysis_path, 'results/models/loocv_baseline_map_IVb_bumps5',comp_add,'.mat'));
likelis_best = loocv_best.likelihoods;

for map_i = 1:length(maps)
    map_index = map_indices(map_i)
    cur_map = maps{map_i};

    %specify nr of bumps
    Nbumps = sum(cur_map(1,:))-1;
    
    disp(['Map - Nbumps ', num2str(Nbumps)])
    
    likelihoods = repmat(-Inf,nsubs,1); %BestLikelihood
    bumps = zeros(pcs,size(cur_map,2)-1,nsubs); %BestBumpMag
    gammas = zeros(size(cur_map,2), 2, size(cur_map,1), nsubs); %BestGamma
    
    parpool(5)
    
    parfor subj = 1:nsubs

        % Select starting parameters for this fold
        new_bumps = optimal_bumps{subj,Nbumps};
        new_gammas = optimal_gammas{subj,Nbumps};
    
        colsums = sum(cur_map,1);
        if ismember(1,colsums)
            have_bump = [];
            ci = 1;
            col = 1;
            % Identify which stages have a different bump per word type
            while ci < length(colsums)
                if colsums(ci) == 1
                    have_bump(end+1) = col;
                    ci = ci + (ncond -1);
                end
                col = col + 1;
                ci = ci + 1;
            end
            
            % Now add for each word type and identified stage extra bump and
            % gamma parameters.
            if ~ismember(1,have_bump)
                tmp_bumps = new_bumps(:,1);
                tmp_gammas = new_gammas(1,:);
            else
                tmp_bumps = repmat(new_bumps(:,1),1,ncond);
                tmp_gammas = repmat(new_gammas(1,:),ncond,1);
            end
            
            for stage_i = 2:(Nbumps+1)
                
                if ~ismember(stage_i,have_bump)
                    if stage_i <= Nbumps % One bump less than number of stages!
                        tmp_bumps = horzcat(tmp_bumps,new_bumps(:,stage_i));
                    end
                    tmp_gammas = vertcat(tmp_gammas,new_gammas(stage_i,:));
                else
                    if stage_i <= Nbumps
                        tmp_bumps = horzcat(tmp_bumps,repmat(new_bumps(:,stage_i),1,ncond));
                    end
                    tmp_gammas = vertcat(tmp_gammas,repmat(new_gammas(stage_i,:),ncond,1));
                end
            end
            new_bumps = tmp_bumps
            new_gammas = tmp_gammas
        end
    
        %select data one
        score_one = score(subjects_varS == subj,:);
    
        %x/y one
        trial_lens_one = trial_lens(subjects_var == subj);
        x_one = [0; cumsum(trial_lens_one)] + 1;
        x_one(end) = [];
        y_one = x_one + trial_lens_one - 1;
        conds_one = conditions(subjects_var == subj);
        
    
        %select data rest
        score_rest = score(subjects_varS ~= subj,:);
    
        %x/y rest
        trial_lens_rest = trial_lens(subjects_var ~= subj); 
        x_rest = [0; cumsum(trial_lens_rest)] + 1;
        x_rest(end) = [];
        y_rest = x_rest + trial_lens_rest - 1;
        conds_rest = conditions(subjects_var ~= subj);
    
        %calculate model based on other subjects
        [~, bumps(:,:,subj), gammas(:,:,:,subj),~]=hsmmBoundCorrsCondsVary(score_rest, conds_rest, ...
                new_bumps,new_gammas, ...
                1,x_rest,y_rest,max_dur,cur_map,c);
            
    
        % fit the model on the left out subject
        [likelihoods(subj), ~, ~, ~] = hsmmBoundCorrsCondsVary(score_one,conds_one, bumps(:,:,subj), gammas(:,:,:,subj),0,x_one,y_one,max_dur, cur_map, c);
    end
    
    % calc final solution based on average
    % subject
    [~,bumps_final,gammas_final, events_final]= hsmmBoundCorrsCondsVary(score,conditions,nanmean(bumps,3),nanmean(gammas,4),1,x,y,max_dur,cur_map,c);
    
    % save to
    save(strcat(analysis_path, 'results/models/loocv_baseline_map_',map_index,'_bumps', num2str(Nbumps),comp_add,'.mat'),...
    'likelihoods', 'bumps', 'gammas', 'bumps_final', 'gammas_final','events_final');
    
    % Check improvement according to current candidate
    improv = sum(likelihoods > likelis_best)
    [p, h] = signtest(likelihoods,likelis_best,'Tail','right')
    
    delete(gcp('nocreate'))


end

% Now test the extra bump maps

maps = {map_VI,map_VII,map_VIII};
map_indices = ["VIc","VIIc","VIIIc"];

for map_i = 1:length(maps)
    map_index = map_indices(map_i)
    cur_map = maps{map_i};

    extra_bump_cond = 1:3;
    extra_bump_cond = extra_bump_cond(cur_map(:,7)' == 1)

    likelihoods = repmat(-Inf,nsubs,1); %BestLikelihood
    bumps = zeros(pcs,size(cur_map,2)-1,nsubs); %BestBumpMag
    gammas = zeros(size(cur_map,2), 2, size(cur_map,1), nsubs); %BestGamma
    
    parpool(5)
    
    parfor subj = 1:nsubs
    
        % Select starting parameters for this fold
        % first pick shared bumps
        new_bumps = optimal_bumps{subj,5};
        new_gammas = optimal_gammas{subj,5};
    
        % Now add initial parameters for the first three stages
        % based on the shared solution
        tmp_bumps = new_bumps(:,1:3);
        tmp_gammas = new_gammas(1:3,:);
    
        % Bump preceding stage 5 (six for two word types but one in this
        % extra bump map) and gammas of stage 4 differ per word types.
        tmp_bumps = horzcat(tmp_bumps,repmat(new_bumps(:,4),1,ncond));
        tmp_gammas = vertcat(tmp_gammas,repmat(new_gammas(4,:),ncond,1));
    
        % Now, we test whether stage 5 exists (only for) a single word type
        % so we insert a copy of the bump preceding stage 5 for that word
        % type only - to see whether it ends up being sufficiently
        % different to justify an extra stage/bump topology. Since the
        % stage falls right between the original stages four/five we use the average of
        % the gammas to initialize the gamma for this extra stage.
        tmp_bumps = horzcat(tmp_bumps,new_bumps(:,4));
        tmp_gammas = vertcat(tmp_gammas,nanmean(new_gammas(4:5,:),1));
    
        % For stage six (five for two word types) we need separate
        % gamma parameters per word type as well as a shared bump preceding
        % stage seven (six for two word types).
        tmp_bumps = horzcat(tmp_bumps,new_bumps(:,5));
        tmp_gammas = vertcat(tmp_gammas,new_gammas(5,:));
    
        % For stage seven (six for two word types) we need a shared
        % gamma parameter.
        tmp_gammas = vertcat(tmp_gammas,new_gammas(6,:));
        new_bumps = tmp_bumps;
        new_gammas = tmp_gammas;
    
        %select data one
        score_one = score(subjects_varS == subj,:);
    
        %x/y one
        trial_lens_one = trial_lens(subjects_var == subj);
        x_one = [0; cumsum(trial_lens_one)] + 1;
        x_one(end) = [];
        y_one = x_one + trial_lens_one - 1;
        conds_one = conditions(subjects_var == subj);
        
        %select data rest
        score_rest = score(subjects_varS ~= subj,:);
    
        %x/y rest
        trial_lens_rest = trial_lens(subjects_var ~= subj); 
        x_rest = [0; cumsum(trial_lens_rest)] + 1;
        x_rest(end) = [];
        y_rest = x_rest + trial_lens_rest - 1;
        conds_rest = conditions(subjects_var ~= subj);
    
        %calculate model based on other subjects
        [~, bumps(:,:,subj), gammas(:,:,:,subj),~]=hsmmBoundCorrsCondsVary(score_rest, conds_rest, ...
                new_bumps,new_gammas, ...
                1,x_rest,y_rest,max_dur,cur_map,c);
            
    
        % fit the model on the left out subject
        [likelihoods(subj), ~, ~, ~] = hsmmBoundCorrsCondsVary(score_one,conds_one, bumps(:,:,subj), gammas(:,:,:,subj),0,x_one,y_one,max_dur, cur_map, c);
    end
    
    % calc final solution based on average
    % subject
    [~,bumps_final,gammas_final, events_final]= hsmmBoundCorrsCondsVary(score,conditions,nanmean(bumps,3),nanmean(gammas,4),1,x,y,max_dur,cur_map,c);
    
    % save to
    save(strcat(analysis_path, 'results/models/loocv_baseline_map_',map_index,'_bumps', num2str(7),comp_add,'.mat'),...
    'likelihoods', 'bumps', 'gammas', 'bumps_final', 'gammas_final','events_final');
    
    % Check improvement according to best map
    improv = sum(likelihoods > likelis_best)
    [p, h] = signtest(likelihoods,likelis_best,'Tail','right')

    delete(gcp('nocreate'))

end

%% Fifth forward pass
% Stage five was split into two for pseudo-words. Hence, we need to
% re-check gammas + bumps.
% We start with separate Gammas:

map_I = [1 1 1 1 0 0 0 1 1;... 
         2 1 1 0 1 0 1 2 1;...
         3 1 1 0 0 1 0 3 1];

map_II = [1 1 1 1 0 0 0 1 1;... 
          1 2 1 0 1 0 1 2 1;...
          1 3 1 0 0 1 0 3 1];

map_III = [1 1 1 1 0 0 0 1 1;... 
           1 1 2 0 1 0 1 2 1;...
           1 1 3 0 0 1 0 3 1];

map_IV = [1 1 1 1 0 0 0 1 1;... 
          1 1 1 0 1 0 1 2 2;...
          1 1 1 0 0 1 0 3 3];

% We also have to check whether bump preceding stage six/seven (i.e.,
% the last bump) should differ per word type:

map_V = [1 1 1 1 0 0 0 1 0 0 1;... 
         1 1 1 0 1 0 1 0 1 0 1;...
         1 1 1 0 0 1 0 0 0 1 1];

maps = {map_I,map_II,map_III,map_IV};
map_indices = ["Id","IId","IIId","IVd"];

% re-load data
load([saveTo strcat('data4HMM_baseline_LD',comp_add,'.mat')]);

%load optimal solutions from fitbumps overall model
res_path = strcat('results/init_params/LD_baseline_all_sols_singleModel/');
load(strcat(analysis_path, res_path, 'optimal_sols_baseline_singleModel',comp_add,'.mat'));

%select data
pcs=10;
c=0;
score = normedscore10(:,1:pcs);

%rebuild x/y
trial_lens = y-x+1; %calc lengths //dur

max_dur = max(trial_lens);
min_dur = min(trial_lens);
max_bumps = floor(min_dur / 5);

% Number of word types
ncond = 3;

% Get likelihood from current best model:
loocv_best = load(strcat(analysis_path, 'results/models/loocv_baseline_map_VIIc_bumps7',comp_add,'.mat'));
likelis_best = loocv_best.likelihoods;

% Start with maps only switching gammas (code remains the same as for
% the second set of maps used in the last pass)

for map_i = 1:length(maps)
    map_index = map_indices(map_i)
    cur_map = maps{map_i};

    likelihoods = repmat(-Inf,nsubs,1); %BestLikelihood
    bumps = zeros(pcs,size(cur_map,2)-1,nsubs); %BestBumpMag
    gammas = zeros(size(cur_map,2), 2, size(cur_map,1), nsubs); %BestGamma
    
    parpool(5)
    
    parfor subj = 1:nsubs
    
        new_bumps = optimal_bumps{subj,5};
        new_gammas = optimal_gammas{subj,5};
        
        % Stages 1-3
        tmp_bumps = new_bumps(:,1:3);
        tmp_gammas = new_gammas(1:3,:);
        
        % Stage four
        tmp_bumps = horzcat(tmp_bumps,repmat(new_bumps(:,4),1,ncond));
        tmp_gammas = vertcat(tmp_gammas,repmat(new_gammas(4,:),ncond,1));
        
        % Stage five (pw only)
        tmp_bumps = horzcat(tmp_bumps,new_bumps(:,4));
        tmp_gammas = vertcat(tmp_gammas,nanmean(new_gammas(4:5,:),1));
        
        % Stage six/five
        tmp_bumps = horzcat(tmp_bumps,new_bumps(:,5));
        tmp_gammas = vertcat(tmp_gammas,new_gammas(5,:));
    
        % Stage seven/six
        tmp_gammas = vertcat(tmp_gammas,new_gammas(6,:));
        new_bumps = tmp_bumps;
        new_gammas = tmp_gammas;
    
        %select data one
        score_one = score(subjects_varS == subj,:);
    
        %x/y one
        trial_lens_one = trial_lens(subjects_var == subj);
        x_one = [0; cumsum(trial_lens_one)] + 1;
        x_one(end) = [];
        y_one = x_one + trial_lens_one - 1;
        conds_one = conditions(subjects_var == subj);
        
        %select data rest
        score_rest = score(subjects_varS ~= subj,:);
    
        %x/y rest
        trial_lens_rest = trial_lens(subjects_var ~= subj); 
        x_rest = [0; cumsum(trial_lens_rest)] + 1;
        x_rest(end) = [];
        y_rest = x_rest + trial_lens_rest - 1;
        conds_rest = conditions(subjects_var ~= subj);
    
        %calculate model based on other subjects
        [~, bumps(:,:,subj), gammas(:,:,:,subj),~]=hsmmBoundCorrsCondsVary(score_rest, conds_rest, ...
                new_bumps,new_gammas, ...
                1,x_rest,y_rest,max_dur,cur_map,c);
            
    
        % fit the model on the left out subject
        [likelihoods(subj), ~, ~, ~] = hsmmBoundCorrsCondsVary(score_one,conds_one, bumps(:,:,subj), gammas(:,:,:,subj),0,x_one,y_one,max_dur, cur_map, c);
    end
    
    % calc final solution based on average
    % subject
    [~,bumps_final,gammas_final, events_final]= hsmmBoundCorrsCondsVary(score,conditions,nanmean(bumps,3),nanmean(gammas,4),1,x,y,max_dur,cur_map,c);
    
    % save to
    save(strcat(analysis_path, 'results/models/loocv_baseline_map_',map_index,'_bumps', num2str(7),comp_add,'.mat'),...
    'likelihoods', 'bumps', 'gammas', 'bumps_final', 'gammas_final','events_final');
    
    % Check improvement according to best map
    improv = sum(likelihoods > likelis_best)
    [p, h] = signtest(likelihoods,likelis_best,'Tail','right')

    delete(gcp('nocreate'))

end

% Now check map V, here the bump preceding the last stage needs to
% differ per word type.

maps = {map_V};
map_indices = ["Vd"];

for map_i = 1:length(maps)
    map_index = map_indices(map_i)
    cur_map = maps{map_i};

    likelihoods = repmat(-Inf,nsubs,1); %BestLikelihood
    bumps = zeros(pcs,size(cur_map,2)-1,nsubs); %BestBumpMag
    gammas = zeros(size(cur_map,2), 2, size(cur_map,1), nsubs); %BestGamma
    
    parpool(5)
    
    parfor subj = 1:nsubs
    
        new_bumps = optimal_bumps{subj,5};
        new_gammas = optimal_gammas{subj,5};
        
        % Stages 1-3
        tmp_bumps = new_bumps(:,1:3);
        tmp_gammas = new_gammas(1:3,:);
        
        % Stage four
        tmp_bumps = horzcat(tmp_bumps,repmat(new_bumps(:,4),1,ncond));
        tmp_gammas = vertcat(tmp_gammas,repmat(new_gammas(4,:),ncond,1));
        
        % Stage five (pw only)
        tmp_bumps = horzcat(tmp_bumps,new_bumps(:,4));
        tmp_gammas = vertcat(tmp_gammas,nanmean(new_gammas(4:5,:),1));
        
        % Stage six/five
        % Here we now have to repeat the last bump 3 times - once
        % per word type. In addition we also need to repeat the gamma
        % parameters - which were different for this stage anyway.
        tmp_bumps = horzcat(tmp_bumps,repmat(new_bumps(:,5),1,ncond));
        tmp_gammas = vertcat(tmp_gammas,repmat(new_gammas(5,:),ncond,1));
    
        % Stage seven/six
        % This only concerns gamma parameters - which remain fixed in this
        % map.
        tmp_gammas = vertcat(tmp_gammas,new_gammas(6,:));
        new_bumps = tmp_bumps;
        new_gammas = tmp_gammas;
    
        %select data one
        score_one = score(subjects_varS == subj,:);
    
        %x/y one
        trial_lens_one = trial_lens(subjects_var == subj);
        x_one = [0; cumsum(trial_lens_one)] + 1;
        x_one(end) = [];
        y_one = x_one + trial_lens_one - 1;
        conds_one = conditions(subjects_var == subj);
        
        %select data rest
        score_rest = score(subjects_varS ~= subj,:);
    
        %x/y rest
        trial_lens_rest = trial_lens(subjects_var ~= subj); 
        x_rest = [0; cumsum(trial_lens_rest)] + 1;
        x_rest(end) = [];
        y_rest = x_rest + trial_lens_rest - 1;
        conds_rest = conditions(subjects_var ~= subj);
    
        %calculate model based on other subjects
        [~, bumps(:,:,subj), gammas(:,:,:,subj),~]=hsmmBoundCorrsCondsVary(score_rest, conds_rest, ...
                new_bumps,new_gammas, ...
                1,x_rest,y_rest,max_dur,cur_map,c);
            
    
        % fit the model on the left out subject
        [likelihoods(subj), ~, ~, ~] = hsmmBoundCorrsCondsVary(score_one,conds_one, bumps(:,:,subj), gammas(:,:,:,subj),0,x_one,y_one,max_dur, cur_map, c);
    end
    
    % calc final solution based on average
    % subject
    [~,bumps_final,gammas_final, events_final]= hsmmBoundCorrsCondsVary(score,conditions,nanmean(bumps,3),nanmean(gammas,4),1,x,y,max_dur,cur_map,c);
    
    % save to
    save(strcat(analysis_path, 'results/models/loocv_baseline_map_',map_index,'_bumps', num2str(7),comp_add,'.mat'),...
    'likelihoods', 'bumps', 'gammas', 'bumps_final', 'gammas_final','events_final');
    
    % Check improvement according to best map
    improv = sum(likelihoods > likelis_best)
    [p, h] = signtest(likelihoods,likelis_best,'Tail','right')

    delete(gcp('nocreate'))

end

%% Sixth Forward pass
% No evidence for more complex model. As a final check we need to test
% whether for the models in which only the gamma parameter varied we also
% need to estimate separate bumps.
map_I = [1 0 0 1 1 1 0 0 0 1 1;... 
         0 1 0 1 1 0 1 0 1 2 1;...
         0 0 1 1 1 0 0 1 0 3 1];

map_II = [1 1 0 0 1 1 0 0 0 1 1;... 
          1 0 1 0 1 0 1 0 1 2 1;...
          1 0 0 1 1 0 0 1 0 3 1];

map_III = [1 1 1 0 0 1 0 0 0 1 1;... 
           1 1 0 1 0 0 1 0 1 2 1;...
           1 1 0 0 1 0 0 1 0 3 1];

maps = {map_I,map_II,map_III};
map_indices = ["Ie","IIe","IIIe"];

% re-load data
load([saveTo strcat('data4HMM_baseline_LD',comp_add,'.mat')]);

%load optimal solutions from fitbumps overall model
res_path = strcat('results/init_params/LD_baseline_all_sols_singleModel/');
load(strcat(analysis_path, res_path, 'optimal_sols_baseline_singleModel',comp_add,'.mat'));

%select data
pcs=10;
c=0;
score = normedscore10(:,1:pcs);

%rebuild x/y
trial_lens = y-x+1; %calc lengths //dur

max_dur = max(trial_lens);
min_dur = min(trial_lens);
max_bumps = floor(min_dur / 5);

% Number of word types
ncond = 3;

% Get likelihood from current best model:
loocv_best = load(strcat(analysis_path, 'results/models/loocv_baseline_map_VIIc_bumps7',comp_add,'.mat'));
likelis_best = loocv_best.likelihoods;

% Code again remains almost the same as for the first set of maps used in
% the last pass). What needs to change is the setup for the first three
% stages - one needs to get separate bumps per word type.

for map_i = 1:length(maps)
    map_index = map_indices(map_i)
    cur_map = maps{map_i};

    likelihoods = repmat(-Inf,nsubs,1); %BestLikelihood
    bumps = zeros(pcs,size(cur_map,2)-1,nsubs); %BestBumpMag
    gammas = zeros(size(cur_map,2), 2, size(cur_map,1), nsubs); %BestGamma
    
    parpool(5)
    
    parfor subj = 1:nsubs
    
        new_bumps = optimal_bumps{subj,5};
        new_gammas = optimal_gammas{subj,5};
        
        % Stages 1-3 - take care of separate bumps per word type
        if map_i == 1
            tmp_bumps = repmat(new_bumps(:,1),1,ncond);
            tmp_gammas = repmat(new_gammas(1,:),ncond,1);
        else
            tmp_bumps = new_bumps(:,1);
            tmp_gammas = new_gammas(1,:);
        end
        
        for stg_i = 2:3
            if map_i == stg_i
                tmp_bumps = horzcat(tmp_bumps,repmat(new_bumps(:,stg_i),1,ncond));
                tmp_gammas = vertcat(tmp_gammas,repmat(new_gammas(stg_i,:),ncond,1));
            else
                tmp_bumps = horzcat(tmp_bumps,new_bumps(:,stg_i));
                tmp_gammas = vertcat(tmp_gammas,new_gammas(stg_i,:));
            end
        end
        
        % Stage four
        tmp_bumps = horzcat(tmp_bumps,repmat(new_bumps(:,4),1,ncond));
        tmp_gammas = vertcat(tmp_gammas,repmat(new_gammas(4,:),ncond,1));
        
        % Stage five (pw only)
        tmp_bumps = horzcat(tmp_bumps,new_bumps(:,4));
        tmp_gammas = vertcat(tmp_gammas,nanmean(new_gammas(4:5,:),1));
        
        % Stage six/five
        tmp_bumps = horzcat(tmp_bumps,new_bumps(:,5));
        tmp_gammas = vertcat(tmp_gammas,new_gammas(5,:));
    
        % Stage seven/six
        tmp_gammas = vertcat(tmp_gammas,new_gammas(6,:));
        new_bumps = tmp_bumps;
        new_gammas = tmp_gammas;
    
        %select data one
        score_one = score(subjects_varS == subj,:);
    
        %x/y one
        trial_lens_one = trial_lens(subjects_var == subj);
        x_one = [0; cumsum(trial_lens_one)] + 1;
        x_one(end) = [];
        y_one = x_one + trial_lens_one - 1;
        conds_one = conditions(subjects_var == subj);
        
        %select data rest
        score_rest = score(subjects_varS ~= subj,:);
    
        %x/y rest
        trial_lens_rest = trial_lens(subjects_var ~= subj); 
        x_rest = [0; cumsum(trial_lens_rest)] + 1;
        x_rest(end) = [];
        y_rest = x_rest + trial_lens_rest - 1;
        conds_rest = conditions(subjects_var ~= subj);
    
        %calculate model based on other subjects
        [~, bumps(:,:,subj), gammas(:,:,:,subj),~]=hsmmBoundCorrsCondsVary(score_rest, conds_rest, ...
                new_bumps,new_gammas, ...
                1,x_rest,y_rest,max_dur,cur_map,c);
            
    
        % fit the model on the left out subject
        [likelihoods(subj), ~, ~, ~] = hsmmBoundCorrsCondsVary(score_one,conds_one, bumps(:,:,subj), gammas(:,:,:,subj),0,x_one,y_one,max_dur, cur_map, c);
    end
    
    % calc final solution based on average
    % subject
    [~,bumps_final,gammas_final, events_final]= hsmmBoundCorrsCondsVary(score,conditions,nanmean(bumps,3),nanmean(gammas,4),1,x,y,max_dur,cur_map,c);
    
    % save to
    save(strcat(analysis_path, 'results/models/loocv_baseline_map_',map_index,'_bumps', num2str(7),comp_add,'.mat'),...
    'likelihoods', 'bumps', 'gammas', 'bumps_final', 'gammas_final','events_final');
    
    % Check improvement according to best map
    improv = sum(likelihoods > likelis_best)
    [p, h] = signtest(likelihoods,likelis_best,'Tail','right')

    delete(gcp('nocreate'))

end

%% Seventh forward pass
% Bump preceding stage four differs between word types. We need to again
% consider new gamma + bump parameters. In addition we need to check
% whether stage four needs to be split in half for random strings and words
% this cannot be the case for pseudo-words since the maximum is seven
% stages.

map_I = [1 1 1 0 0 1 0 0 0 1 1;... 
         2 1 0 1 0 0 1 0 1 2 1;...
         3 1 0 0 1 0 0 1 0 3 1];

map_II = [1 1 1 0 0 1 0 0 0 1 1;... 
          1 2 0 1 0 0 1 0 1 2 1;...
          1 3 0 0 1 0 0 1 0 3 1];

map_III = [1 1 1 0 0 1 0 0 0 1 1;... 
           1 1 0 1 0 0 1 0 1 2 2;...
           1 1 0 0 1 0 0 1 0 3 3];

% Again need to test whether last bump should differ between word types
map_IV = [1 1 1 0 0 1 0 0 0 1 0 0 1;... 
          1 1 0 1 0 0 1 0 1 0 1 0 1;...
          1 1 0 0 1 0 0 1 0 0 0 1 1];

% Finally - we need to test whether stage four should be split for words
% or random strings.

map_V = [1 1 1 0 0 1 1 0 0 0 1 1;... % Extra stage for words
         1 1 0 1 0 0 0 1 0 1 2 1;...
         1 1 0 0 1 0 0 0 1 0 3 1];

map_VI = [1 1 1 0 0 0 1 0 0 0 1 1;... % Extra stage for random strings
          1 1 0 1 0 0 0 1 0 1 2 1;...
          1 1 0 0 1 1 0 0 1 0 3 1];


% We start with testing the Gamma maps again
maps = {map_I,map_II,map_III};
map_indices = ["If","IIf","IIIf"];

% re-load data
load([saveTo strcat('data4HMM_baseline_LD',comp_add,'.mat')]);

%load optimal solutions from fitbumps overall model
res_path = strcat('results/init_params/LD_baseline_all_sols_singleModel/');
load(strcat(analysis_path, res_path, 'optimal_sols_baseline_singleModel',comp_add,'.mat'));

%select data
pcs=10;
c=0;
score = normedscore10(:,1:pcs);

%rebuild x/y
trial_lens = y-x+1; %calc lengths //dur

max_dur = max(trial_lens);
min_dur = min(trial_lens);
max_bumps = floor(min_dur / 5);

% Number of word types
ncond = 3;

% Get likelihood from current best model:
loocv_best = load(strcat(analysis_path, 'results/models/loocv_baseline_map_IIIe_bumps7',comp_add,'.mat'));
likelis_best = loocv_best.likelihoods;

% Code again remains almost the same as for the last set.
for map_i = 1:length(maps)
    map_index = map_indices(map_i)
    cur_map = maps{map_i};

    likelihoods = repmat(-Inf,nsubs,1); %BestLikelihood
    bumps = zeros(pcs,size(cur_map,2)-1,nsubs); %BestBumpMag
    gammas = zeros(size(cur_map,2), 2, size(cur_map,1), nsubs); %BestGamma
    
    parpool(5)
    
    parfor subj = 1:nsubs
    
        new_bumps = optimal_bumps{subj,5};
        new_gammas = optimal_gammas{subj,5};
        
        % Stages 1-2
        tmp_bumps = new_bumps(:,1:2);
        tmp_gammas = new_gammas(1:2,:);
        
        % Stage three - separate bumps preceding stage four per word type
        tmp_bumps = horzcat(tmp_bumps,repmat(new_bumps(:,3),1,ncond));
        tmp_gammas = vertcat(tmp_gammas,repmat(new_gammas(3,:),ncond,1));
            
        % Stage four
        tmp_bumps = horzcat(tmp_bumps,repmat(new_bumps(:,4),1,ncond));
        tmp_gammas = vertcat(tmp_gammas,repmat(new_gammas(4,:),ncond,1));
        
        % Stage five (pw only)
        tmp_bumps = horzcat(tmp_bumps,new_bumps(:,4));
        tmp_gammas = vertcat(tmp_gammas,nanmean(new_gammas(4:5,:),1));
        
        % Stage six/five
        tmp_bumps = horzcat(tmp_bumps,new_bumps(:,5));
        tmp_gammas = vertcat(tmp_gammas,new_gammas(5,:));
    
        % Stage seven/six
        tmp_gammas = vertcat(tmp_gammas,new_gammas(6,:));
        new_bumps = tmp_bumps;
        new_gammas = tmp_gammas;
    
        %select data one
        score_one = score(subjects_varS == subj,:);
    
        %x/y one
        trial_lens_one = trial_lens(subjects_var == subj);
        x_one = [0; cumsum(trial_lens_one)] + 1;
        x_one(end) = [];
        y_one = x_one + trial_lens_one - 1;
        conds_one = conditions(subjects_var == subj);
        
        %select data rest
        score_rest = score(subjects_varS ~= subj,:);
    
        %x/y rest
        trial_lens_rest = trial_lens(subjects_var ~= subj); 
        x_rest = [0; cumsum(trial_lens_rest)] + 1;
        x_rest(end) = [];
        y_rest = x_rest + trial_lens_rest - 1;
        conds_rest = conditions(subjects_var ~= subj);
    
        %calculate model based on other subjects
        [~, bumps(:,:,subj), gammas(:,:,:,subj),~]=hsmmBoundCorrsCondsVary(score_rest, conds_rest, ...
                new_bumps,new_gammas, ...
                1,x_rest,y_rest,max_dur,cur_map,c);
            
    
        % fit the model on the left out subject
        [likelihoods(subj), ~, ~, ~] = hsmmBoundCorrsCondsVary(score_one,conds_one, bumps(:,:,subj), gammas(:,:,:,subj),0,x_one,y_one,max_dur, cur_map, c);
    end
    
    % calc final solution based on average
    % subject
    [~,bumps_final,gammas_final, events_final]= hsmmBoundCorrsCondsVary(score,conditions,nanmean(bumps,3),nanmean(gammas,4),1,x,y,max_dur,cur_map,c);
    
    % save to
    save(strcat(analysis_path, 'results/models/loocv_baseline_map_',map_index,'_bumps', num2str(7),comp_add,'.mat'),...
    'likelihoods', 'bumps', 'gammas', 'bumps_final', 'gammas_final','events_final');
    
    % Check improvement according to best map
    improv = sum(likelihoods > likelis_best)
    [p, h] = signtest(likelihoods,likelis_best,'Tail','right')

    delete(gcp('nocreate'))

end

% Now maps IV - with separate last bumps
maps = {map_IV};
map_indices = ["IVf"];

for map_i = 1:length(maps)
    map_index = map_indices(map_i)
    cur_map = maps{map_i};

    likelihoods = repmat(-Inf,nsubs,1); %BestLikelihood
    bumps = zeros(pcs,size(cur_map,2)-1,nsubs); %BestBumpMag
    gammas = zeros(size(cur_map,2), 2, size(cur_map,1), nsubs); %BestGamma
    
    parpool(5)
    
    parfor subj = 1:nsubs
    
        new_bumps = optimal_bumps{subj,5};
        new_gammas = optimal_gammas{subj,5};
        
        % Stages 1-2
        tmp_bumps = new_bumps(:,1:2);
        tmp_gammas = new_gammas(1:2,:);
        
        % Stage three - separate bumps preceding stage four per word type
        tmp_bumps = horzcat(tmp_bumps,repmat(new_bumps(:,3),1,ncond));
        tmp_gammas = vertcat(tmp_gammas,repmat(new_gammas(3,:),ncond,1));
            
        % Stage four
        tmp_bumps = horzcat(tmp_bumps,repmat(new_bumps(:,4),1,ncond));
        tmp_gammas = vertcat(tmp_gammas,repmat(new_gammas(4,:),ncond,1));
        
        % Stage five (pw only)
        tmp_bumps = horzcat(tmp_bumps,new_bumps(:,4));
        tmp_gammas = vertcat(tmp_gammas,nanmean(new_gammas(4:5,:),1));
        
        % Stage six/five
        % Here we now have to repeat the last bump 3 times - once
        % per word type. In addition we also need to repeat the gamma
        % parameters - which were different for this stage anyway.
        tmp_bumps = horzcat(tmp_bumps,repmat(new_bumps(:,5),1,ncond));
        tmp_gammas = vertcat(tmp_gammas,repmat(new_gammas(5,:),ncond,1));
    
        % Stage seven/six
        % This only concerns gamma parameters - which remain fixed in this
        % map.
        tmp_gammas = vertcat(tmp_gammas,new_gammas(6,:));
        new_bumps = tmp_bumps;
        new_gammas = tmp_gammas;
    
        %select data one
        score_one = score(subjects_varS == subj,:);
    
        %x/y one
        trial_lens_one = trial_lens(subjects_var == subj);
        x_one = [0; cumsum(trial_lens_one)] + 1;
        x_one(end) = [];
        y_one = x_one + trial_lens_one - 1;
        conds_one = conditions(subjects_var == subj);
        
        %select data rest
        score_rest = score(subjects_varS ~= subj,:);
    
        %x/y rest
        trial_lens_rest = trial_lens(subjects_var ~= subj); 
        x_rest = [0; cumsum(trial_lens_rest)] + 1;
        x_rest(end) = [];
        y_rest = x_rest + trial_lens_rest - 1;
        conds_rest = conditions(subjects_var ~= subj);
    
        %calculate model based on other subjects
        [~, bumps(:,:,subj), gammas(:,:,:,subj),~]=hsmmBoundCorrsCondsVary(score_rest, conds_rest, ...
                new_bumps,new_gammas, ...
                1,x_rest,y_rest,max_dur,cur_map,c);
            
    
        % fit the model on the left out subject
        [likelihoods(subj), ~, ~, ~] = hsmmBoundCorrsCondsVary(score_one,conds_one, bumps(:,:,subj), gammas(:,:,:,subj),0,x_one,y_one,max_dur, cur_map, c);
    end
    
    % calc final solution based on average
    % subject
    [~,bumps_final,gammas_final, events_final]= hsmmBoundCorrsCondsVary(score,conditions,nanmean(bumps,3),nanmean(gammas,4),1,x,y,max_dur,cur_map,c);
    
    % save to
    save(strcat(analysis_path, 'results/models/loocv_baseline_map_',map_index,'_bumps', num2str(7),comp_add,'.mat'),...
    'likelihoods', 'bumps', 'gammas', 'bumps_final', 'gammas_final','events_final');
    
    % Check improvement according to best map
    improv = sum(likelihoods > likelis_best)
    [p, h] = signtest(likelihoods,likelis_best,'Tail','right')

    delete(gcp('nocreate'))

end

% Finally, the two maps with extra bumps for words and random strings
maps = {map_V,map_VI};
map_indices = ["Vf","VIf"];

for map_i = 1:length(maps)
    map_index = map_indices(map_i)
    cur_map = maps{map_i};

    likelihoods = repmat(-Inf,nsubs,1); %BestLikelihood
    bumps = zeros(pcs,size(cur_map,2)-1,nsubs); %BestBumpMag
    gammas = zeros(size(cur_map,2), 2, size(cur_map,1), nsubs); %BestGamma
    
    parpool(5)
    
    parfor subj = 1:nsubs
    
        new_bumps = optimal_bumps{subj,5};
        new_gammas = optimal_gammas{subj,5};
        
        % Stages 1-2
        tmp_bumps = new_bumps(:,1:2);
        tmp_gammas = new_gammas(1:2,:);
        
        % Stage three - separate bumps preceding stage four per word type
        tmp_bumps = horzcat(tmp_bumps,repmat(new_bumps(:,3),1,ncond));
        tmp_gammas = vertcat(tmp_gammas,repmat(new_gammas(3,:),ncond,1));

        % Stage four (words or random strings only).
        % we need to test whether stage four can be split in half for words
        % or random strings. So we copy bump 3 (preceding stage four) for
        % the word type to be tested and again initialize the gamma for
        % this extra stage based on the average of what was stage 3 and 4.
        tmp_bumps = horzcat(tmp_bumps,new_bumps(:,3));
        tmp_gammas = vertcat(tmp_gammas,nanmean(new_gammas(3:4,:),1));

        % Stage four (five for words or random strings)
        tmp_bumps = horzcat(tmp_bumps,repmat(new_bumps(:,4),1,ncond));
        tmp_gammas = vertcat(tmp_gammas,repmat(new_gammas(4,:),ncond,1));
        
        % Stage five (pw only)
        tmp_bumps = horzcat(tmp_bumps,new_bumps(:,4));
        tmp_gammas = vertcat(tmp_gammas,nanmean(new_gammas(4:5,:),1));
        
        % Stage six/five
        tmp_bumps = horzcat(tmp_bumps,new_bumps(:,5));
        tmp_gammas = vertcat(tmp_gammas,new_gammas(5,:));
    
        % Stage seven/six
        tmp_gammas = vertcat(tmp_gammas,new_gammas(6,:));
        new_bumps = tmp_bumps;
        new_gammas = tmp_gammas;
    
        %select data one
        score_one = score(subjects_varS == subj,:);
    
        %x/y one
        trial_lens_one = trial_lens(subjects_var == subj);
        x_one = [0; cumsum(trial_lens_one)] + 1;
        x_one(end) = [];
        y_one = x_one + trial_lens_one - 1;
        conds_one = conditions(subjects_var == subj);
        
        %select data rest
        score_rest = score(subjects_varS ~= subj,:);
    
        %x/y rest
        trial_lens_rest = trial_lens(subjects_var ~= subj); 
        x_rest = [0; cumsum(trial_lens_rest)] + 1;
        x_rest(end) = [];
        y_rest = x_rest + trial_lens_rest - 1;
        conds_rest = conditions(subjects_var ~= subj);
    
        %calculate model based on other subjects
        [~, bumps(:,:,subj), gammas(:,:,:,subj),~]=hsmmBoundCorrsCondsVary(score_rest, conds_rest, ...
                new_bumps,new_gammas, ...
                1,x_rest,y_rest,max_dur,cur_map,c);
            
    
        % fit the model on the left out subject
        [likelihoods(subj), ~, ~, ~] = hsmmBoundCorrsCondsVary(score_one,conds_one, bumps(:,:,subj), gammas(:,:,:,subj),0,x_one,y_one,max_dur, cur_map, c);
    end
    
    % calc final solution based on average
    % subject
    [~,bumps_final,gammas_final, events_final]= hsmmBoundCorrsCondsVary(score,conditions,nanmean(bumps,3),nanmean(gammas,4),1,x,y,max_dur,cur_map,c);
    
    % save to
    save(strcat(analysis_path, 'results/models/loocv_baseline_map_',map_index,'_bumps', num2str(7),comp_add,'.mat'),...
    'likelihoods', 'bumps', 'gammas', 'bumps_final', 'gammas_final','events_final');
    
    % Check improvement according to best map
    improv = sum(likelihoods > likelis_best)
    [p, h] = signtest(likelihoods,likelis_best,'Tail','right')

    delete(gcp('nocreate'))

end

%% Eight forward pass
% No evidence for more complex model - extra stages for words and random
% strings were not supported by the data. We now need to test extra bumps
% for the previous gamma only models again.

map_I = [1 0 0 1 1 0 0 1 0 0 0 1 1;... 
         0 1 0 1 0 1 0 0 1 0 1 2 1;...
         0 0 1 1 0 0 1 0 0 1 0 3 1];

map_II = [1 1 0 0 1 0 0 1 0 0 0 1 1;... 
          1 0 1 0 0 1 0 0 1 0 1 2 1;...
          1 0 0 1 0 0 1 0 0 1 0 3 1];

maps = {map_I,map_II};
map_indices = ["Ig","IIg"];

% re-load data
load([saveTo strcat('data4HMM_baseline_LD',comp_add,'.mat')]);

%load optimal solutions from fitbumps overall model
res_path = strcat('results/init_params/LD_baseline_all_sols_singleModel/');
load(strcat(analysis_path, res_path, 'optimal_sols_baseline_singleModel',comp_add,'.mat'));

%select data
pcs=10;
c=0;
score = normedscore10(:,1:pcs);

%rebuild x/y
trial_lens = y-x+1; %calc lengths //dur

max_dur = max(trial_lens);
min_dur = min(trial_lens);
max_bumps = floor(min_dur / 5);

% Number of word types
ncond = 3;

% Get likelihood from current best model:
loocv_best = load(strcat(analysis_path, 'results/models/loocv_baseline_map_IIIe_bumps7',comp_add,'.mat'));
likelis_best = loocv_best.likelihoods;

% Code again remains almost the same as the first comparison code
% from the last set.
for map_i = 1:length(maps)
    map_index = map_indices(map_i)
    cur_map = maps{map_i};

    likelihoods = repmat(-Inf,nsubs,1); %BestLikelihood
    bumps = zeros(pcs,size(cur_map,2)-1,nsubs); %BestBumpMag
    gammas = zeros(size(cur_map,2), 2, size(cur_map,1), nsubs); %BestGamma
    
    parpool(5)
    
    parfor subj = 1:nsubs
    
        new_bumps = optimal_bumps{subj,5};
        new_gammas = optimal_gammas{subj,5};
        
        % Stages 1-2 - take care of separate bumps per word type
        if map_i == 1
            tmp_bumps = repmat(new_bumps(:,1),1,ncond);
            tmp_gammas = repmat(new_gammas(1,:),ncond,1);

            % Stage 2
            tmp_bumps = horzcat(tmp_bumps,new_bumps(:,2));
            tmp_gammas = vertcat(tmp_gammas,new_gammas(2,:));
        else
            tmp_bumps = new_bumps(:,1);
            tmp_gammas = new_gammas(1,:);
            
            % Stage 2
            tmp_bumps = horzcat(tmp_bumps,repmat(new_bumps(:,2),1,ncond));
            tmp_gammas = vertcat(tmp_gammas,repmat(new_gammas(2,:),ncond,1));
        end
  
        % Stage three - separate bumps preceding stage four per word type
        tmp_bumps = horzcat(tmp_bumps,repmat(new_bumps(:,3),1,ncond));
        tmp_gammas = vertcat(tmp_gammas,repmat(new_gammas(3,:),ncond,1));
            
        % Stage four
        tmp_bumps = horzcat(tmp_bumps,repmat(new_bumps(:,4),1,ncond));
        tmp_gammas = vertcat(tmp_gammas,repmat(new_gammas(4,:),ncond,1));
        
        % Stage five (pw only)
        tmp_bumps = horzcat(tmp_bumps,new_bumps(:,4));
        tmp_gammas = vertcat(tmp_gammas,nanmean(new_gammas(4:5,:),1));
        
        % Stage six/five
        tmp_bumps = horzcat(tmp_bumps,new_bumps(:,5));
        tmp_gammas = vertcat(tmp_gammas,new_gammas(5,:));
    
        % Stage seven/six
        tmp_gammas = vertcat(tmp_gammas,new_gammas(6,:));
        new_bumps = tmp_bumps;
        new_gammas = tmp_gammas;
    
        %select data one
        score_one = score(subjects_varS == subj,:);
    
        %x/y one
        trial_lens_one = trial_lens(subjects_var == subj);
        x_one = [0; cumsum(trial_lens_one)] + 1;
        x_one(end) = [];
        y_one = x_one + trial_lens_one - 1;
        conds_one = conditions(subjects_var == subj);
        
        %select data rest
        score_rest = score(subjects_varS ~= subj,:);
    
        %x/y rest
        trial_lens_rest = trial_lens(subjects_var ~= subj); 
        x_rest = [0; cumsum(trial_lens_rest)] + 1;
        x_rest(end) = [];
        y_rest = x_rest + trial_lens_rest - 1;
        conds_rest = conditions(subjects_var ~= subj);
    
        %calculate model based on other subjects
        [~, bumps(:,:,subj), gammas(:,:,:,subj),~]=hsmmBoundCorrsCondsVary(score_rest, conds_rest, ...
                new_bumps,new_gammas, ...
                1,x_rest,y_rest,max_dur,cur_map,c);
            
    
        % fit the model on the left out subject
        [likelihoods(subj), ~, ~, ~] = hsmmBoundCorrsCondsVary(score_one,conds_one, bumps(:,:,subj), gammas(:,:,:,subj),0,x_one,y_one,max_dur, cur_map, c);
    end
    
    % calc final solution based on average
    % subject
    [~,bumps_final,gammas_final, events_final]= hsmmBoundCorrsCondsVary(score,conditions,nanmean(bumps,3),nanmean(gammas,4),1,x,y,max_dur,cur_map,c);
    
    % save to
    save(strcat(analysis_path, 'results/models/loocv_baseline_map_',map_index,'_bumps', num2str(7),comp_add,'.mat'),...
    'likelihoods', 'bumps', 'gammas', 'bumps_final', 'gammas_final','events_final');
    
    % Check improvement according to best map
    improv = sum(likelihoods > likelis_best)
    [p, h] = signtest(likelihoods,likelis_best,'Tail','right')

    delete(gcp('nocreate'))

end
%% Export trial-level estimates from final model
% No evidence for more complex model, so we can export the current one!

% re-load data
load([saveTo strcat('data4HMM_baseline_LD',comp_add,'.mat')]);

% Select final map
map_index = "IIIe";

pcs=10;
c=0;

%select data
score = normedscore10(:,1:pcs);

%rebuild x/y
trial_lens = y-x+1; %calc lengths

max_dur = max(trial_lens);   % duration of each trial
min_dur = min(trial_lens);
max_bumps = floor(min_dur / 5);

Nbumps=7;
ncond=3;

%only get the mapped solution
loocv = load(strcat(analysis_path, 'results/models/loocv_baseline_map_',map_index,'_bumps', num2str(Nbumps),comp_add,'.mat'));
likelis = loocv.likelihoods;

tmp = loocv;

%Plot recovered bump topologies
load(strcat(saveTo,'chanlocs.mat'));
best_bumps = tmp.bumps_final;

topo = reconstruct(best_bumps,coeff10,latent10,mean(data));  % topologies

Nbumps = size(best_bumps,2);
figure('position', [10 10 240*Nbumps 200]);
for i = 1:Nbumps
    subplot(1,Nbumps,i)
    topoplot(topo(i,:), chanlocs, 'maplimits', [-12 12],'style', 'map', 'electrodes','off', 'shading', 'interp','whitebk','on');
    
    if(i == Nbumps)
        cbar;
    end
end

set(gcf,'Color',[1 1 1]);
set(gcf, 'PaperPositionMode', 'auto');
print('-dpng', strcat(analysis_path, plot_path, 'topo_map_',map_index,'_bumps_', num2str(Nbumps),comp_add,'.png'));
close

%Write trial-level results

%calc durations based on event probs

%first, calc eventProbs
best_events = tmp.events_final; %samples x trials x params

nTrials = size(best_events,2);
maxDur = size(best_events,1);

%first get onsets + response
Nbumps = 6;
onsets = zeros(nTrials, Nbumps); %trials x bumps

for tr = 1:nTrials
    tmp = (1:maxDur) * squeeze(best_events(:,tr,:));
    if sum(tmp>0) > 5
        onsets(tr,:) = tmp(tmp>0); %only take existing bumps
    else
        onsets(tr,[1,2,3,4,6]) = tmp(tmp>0); %only take existing bumps
    end
end
onsets(:,end+1) = y-x; % add end trials
onsets = onsets .* 10; % got to ms

durations = zeros(nTrials, 7);
durations(:,1) = onsets(:,1);

for tr = 1:nTrials
    tr_onsets = onsets(tr,onsets(tr,:) > 0);
    tr_onset_diff = tr_onsets(2:end) - tr_onsets(1:end-1);

    if length(tr_onsets) > 6
        durations(tr,2:end) = tr_onset_diff;
    else
        durations(tr,[2,3,4,6,7]) = tr_onset_diff;
    end
end

% export durations and onset for trial-level analysis
trial_dur_exp = table(durations,onsets,subjects_var,trials_var,rts_var,blocks_var,conditions,stims_var);
writetable(trial_dur_exp,strcat(analysis_path,'results/trial_dat_hsmm_final',comp_add,'.csv'));