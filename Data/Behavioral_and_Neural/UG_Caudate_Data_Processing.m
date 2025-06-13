%%UG behavior and neural data.

clear all, close all, clc
%% Set File Paths

%%upload raw data
UG_Caudate_Raw_Data_path = '';
fpath = UG_Caudate_Raw_Data_path;
load(fpath);

%%%% the function shadedErrorBar is used for plotting timeseries neurotransmitter
%%%% estimates and must be added for plots to work. https://github.com/raacampbell/shadedErrorBar
shadedErrorBar_path = '';
addpath(shadedErrorBar_path)

%%%%save path for table to be used in R studio
save_path ='';
save_filename = 'UG_Caudate_Table_R'; 
full_path = fullfile(save_path, sprintf('%s.csv', save_filename));  % Combine file path and file name
%% Windows for smoothing, zscoring, AUC
task_event = 8:10; %%%Time series DA, NA and 5HT during offer presentation
winS = 10; %%%smoothing window 1 second
z_baseline = 1:71; %%%zscore baseline window 7 seconds
analyze_window = 32:71; %%%AUC window 4 seconds

%% Detrending

for i = 1:length(patients_table_behavior(:,1))
    for j = task_event
        for z = 1:length(patients_table_behavior{i, j}(:,1))
        if any(isnan(patients_table_behavior{i, j}(z,:)))
            continue;
        else
            patients_table_behavior{i, j}(z,:) = detrend(patients_table_behavior{i, j}(z,:), 'linear');

        end
        end
    end
end
%% Smoothing
for i = 1:length(patients_table_behavior(:,1))
    for j = task_event
        for c = 1:30
            for z = winS:71
        patients_table_behavior{i, j}(c, z) = mean(patients_table_behavior{i, j}(c, (z+1-winS):z), 2, 'omitnan');
        patients_table_behavior{i, j}(c,1:winS) = NaN;
            end
        end
    end
end
%% ZScore
for i = 1:length(patients_table_behavior(:,1))
    for j = task_event
        for c = 1:30
        patients_table_behavior_mu{i, j}(c, 1) = mean(patients_table_behavior{i, j}(c,z_baseline), 2, "omitnan");
        patients_table_behavior_sd{i, j}(c, 1) = std(patients_table_behavior{i, j}(c,z_baseline), 0, 2, "omitnan");
        patients_table_behavior{i, j}(c,:) = (patients_table_behavior{i, j}(c,:) - patients_table_behavior_mu{i, j}(c))/patients_table_behavior_sd{i, j}(c);
        end
    end
end

%% 

patients_behavior_PD = patients_table_behavior(strcmp(patients_table_behavior(:, 2), 'PD'), :);

patients_behavior_ET = patients_table_behavior(strcmp(patients_table_behavior(:, 2), 'ET'), :);


%% Bayesian Ideal Observer Model

mean_init  = 10; 
k_init     = 4; 
v_init     = 10; 
sigma_init = 4; 

high_mean = 10; 
low_mean  = 4; 

nParticipants = length(patients_table_behavior(:,1));
nTrials       = 30;
proposal_mat  = zeros(nTrials,nParticipants);

norm    = zeros(nParticipants,nTrials);
norm_pe = zeros(nParticipants,nTrials);
var     = zeros(nParticipants,nTrials);
var_pe  = zeros(nParticipants,nTrials);

pos_normpe = struct([]);
neg_normpe = struct([]);
pos_varpe  = struct([]);
neg_varpe  = struct([]);

for p = 1:nParticipants
    for t = 1:nTrials
        proposal_mat(t,p) = patients_table_behavior{p,5}(t,1);
    end

    kt = zeros(1,nTrials);
    vt = zeros(1,nTrials);

    for t = 1:nTrials
        if t == 1
            kt(t) = k_init + 1;
            vt(t) = v_init + 1;
            norm(p,t)    = ((k_init/kt(t)) * mean_init) + ((1/kt(t)) * proposal_mat(t,p));
            norm_pe(p,t) = proposal_mat(t,p) - mean_init;
            var(p,t)     = ((v_init*sigma_init) + ( (k_init/kt(t)) * (proposal_mat(t,p) - mean_init)^2) ) / vt(t);
            var_pe(p,t)  = norm_pe(p,t)^2;
        else
            kt(t) = kt(t-1) + 1;
            vt(t) = vt(t-1) + 1;
            norm(p,t)    = ((kt(t-1)/kt(t)) * norm(p,t-1)) + ((1/kt(t)) * proposal_mat(t,p));
            norm_pe(p,t) = proposal_mat(t,p) - norm(p,t-1);
            var(p,t)     = ( (vt(t-1) * var(p,t-1)) + ( (kt(t-1)/kt(t)) * (proposal_mat(t,p) - norm(p,t-1))^2) ) / vt(t);
            if t == 2
                var_pe(p,t) = norm_pe(p,t)^2 - (vt(t-1)/v_init) * var(p,t-1);
            else
                var_pe(p,t) = norm_pe(p,t)^2 - (vt(t-1)/vt(t-2)) * var(p,t-1);
            end
        end
        if norm_pe(p,t) > 0
            pos_normpe{p,1}(1,t) = t;
        else
            neg_normpe{p,1}(1,t) = t;
        end
        if var_pe(p,t) > 0
            pos_varpe{p,1}(1,t) = t;
        else
            neg_varpe{p,1}(1,t) = t;
        end
    end
    pos_normpe{p,1} = pos_normpe{p,1}(pos_normpe{p,1} ~= 0);
    neg_normpe{p,1} = neg_normpe{p,1}(neg_normpe{p,1} ~= 0);
    pos_varpe{p,1}  = pos_varpe{p,1}(pos_varpe{p,1} ~= 0);
    neg_varpe{p,1}  = neg_varpe{p,1}(neg_varpe{p,1} ~= 0);
end

norm_ET_pe = norm_pe(1:6,:);
norm_PD_pe = norm_pe(7:19,:);

%% ET Offer High Distribution
b1 = [1:15]; %high block 1
b2 = [16:30];%high block 2

HL_ET_idx = find(strcmp(patients_behavior_ET(:,4), "HL"));
LH_ET_idx = find(strcmp(patients_behavior_ET(:,4), "LH"));

for i = 1:length(HL_ET_idx)
    idx = HL_ET_idx(i);  
    ET_caudate_high1_offers(:,i) = patients_behavior_ET{idx,5}(b1,1);
end

for i = 1:length(LH_ET_idx)
    idx = LH_ET_idx(i);  
    ET_caudate_high2_offers(:,i) = patients_behavior_ET{idx,5}(b2,1);
end

ET_caudate_high_offers= [ET_caudate_high1_offers, ET_caudate_high2_offers] ;
ET_caudate_high_offers = ET_caudate_high_offers(:);

ET_caudate_high_offers_mean = mean(ET_caudate_high_offers, 'all', 'omitnan');
ET_caudate_high_offers_sd = std(ET_caudate_high_offers, 0, 'all', 'omitnan');
%% ET Offer Low Distribution

for i = 1:length(HL_ET_idx)
    idx = HL_ET_idx(i);  
    ET_caudate_low2_offers(:,i) = patients_behavior_ET{idx,5}(b2,1);
end

for i = 1:length(LH_ET_idx)
    idx = LH_ET_idx(i);  
    ET_caudate_low1_offers(:,i) = patients_behavior_ET{idx,5}(b1,1);
end

ET_caudate_low_offers= [ET_caudate_low2_offers, ET_caudate_low1_offers] ;
ET_caudate_low_offers = ET_caudate_low_offers(:);

ET_caudate_low_offers_mean = mean(ET_caudate_low_offers, 'all', 'omitnan');
ET_caudate_low_offers_sd = std(ET_caudate_low_offers, 0, 'all', 'omitnan');


%% PD Offer Low Distribution

HL_PD_idx = find(strcmp(patients_behavior_PD(:,4), "HL"));
LH_PD_idx = find(strcmp(patients_behavior_PD(:,4), "LH"));


for i = 1:length(HL_PD_idx)
    idx = HL_PD_idx(i);  % Actual row in the original table
    PD_caudate_high1_offers(:,i) = patients_behavior_PD{idx,5}(b1,1);
end

% Now loop over just those rows
for i = 1:length(LH_PD_idx)
    idx = LH_PD_idx(i);  % Actual row in the original table
    PD_caudate_high2_offers(:,i) = patients_behavior_PD{idx,5}(b2,1);
end

PD_caudate_high_offers= [PD_caudate_high1_offers, PD_caudate_high2_offers] ;
PD_caudate_high_offers = PD_caudate_high_offers(:);

PD_caudate_high_offers_mean = mean(PD_caudate_high_offers, 'all', 'omitnan');
PD_caudate_high_offers_sd = std(PD_caudate_high_offers, 0, 'all', 'omitnan');
%% PD Offer Low Distribution

for i = 1:length(HL_PD_idx)
    idx = HL_PD_idx(i);  % Actual row in the original table
    PD_caudate_low2_offers(:,i) = patients_behavior_PD{idx,5}(b2,1);
end

% Now loop over just those rows
for i = 1:length(LH_PD_idx)
    idx = LH_PD_idx(i);  % Actual row in the original table
    PD_caudate_low1_offers(:,i) = patients_behavior_PD{idx,5}(b1,1);
end

PD_caudate_low_offers= [PD_caudate_low2_offers, PD_caudate_low1_offers] ;
PD_caudate_low_offers = PD_caudate_low_offers(:);

PD_caudate_low_offers_mean = mean(PD_caudate_low_offers, 'all', 'omitnan');
PD_caudate_low_offers_sd = std(PD_caudate_low_offers, 0, 'all', 'omitnan');


%% Emotional Ratings to Norm Prediction Errors

norm_pe_with_emotion = cell(length(patients_table_behavior(:,1)), 1);

for i = 1:length(patients_table_behavior(:,1))
    % Extract norm_pe values for the i-th subject
    subject_norm_pe = norm_pe(i, :);  % 1x30 array, one value for each trial

    % Extract emotional ratings from patients_table_behavior{i, 6}
    subject_emotional_ratings = patients_table_behavior{i, 6}(:, 1);  % 30x1 array

    % Combine norm_pe values and emotional ratings for each trial
    combined_trials = [subject_norm_pe(:), subject_emotional_ratings];  % 30x2 matrix

    norm_pe_with_emotion{i} = combined_trials;
end


bar_start = -5;
num_ratings = 9;
average_norm_pe = nan(1, num_ratings);
sem_norm_pe = nan(1, num_ratings);

% Loop through each emotional rating (1 to 9)
for rating = 1:num_ratings
    norm_pe_values_for_rating = [];

    for i = 1:length(patients_table_behavior(:,1))  % Loop through each subject
        combined_trials = norm_pe_with_emotion{i};

        matching_indices = combined_trials(:, 2) == rating;
        norm_pe_values_for_rating = [norm_pe_values_for_rating; combined_trials(matching_indices, 1)];
    end

    average_norm_pe(rating) = mean(norm_pe_values_for_rating, 'omitnan');
    sem_norm_pe(rating) = std(norm_pe_values_for_rating, 'omitnan') / sqrt(sum(~isnan(norm_pe_values_for_rating)));
end

%% %% PD
pos_norm_PDpe_NTs = cell(size(patients_behavior_PD, 1), 3);
neg_norm_PDpe_NTs = cell(size(patients_behavior_PD, 1), 3);
pos_norm_PDpe = pos_normpe(7:19,1);
neg_norm_PDpe = neg_normpe(7:19,1);

for i = 1:size(patients_behavior_PD, 1)
    % Get the index values for pos_norm_PDpe and neg_norm_PDpe for this subject
    pos_indexes = pos_norm_PDpe{i};
    neg_indexes = neg_norm_PDpe{i};
    
    for col = task_event
        data = patients_behavior_PD{i, col}(:, 1:71);
        
        col_new = col - 7;  % This will map columns 8:9 which is offer reveal DA, NE and 5HT to -> 1:3
        
        % Extract trials based on pos_norm_PDpe indexes and store in corresponding column
        if ~isempty(pos_indexes)
            pos_norm_PDpe_NTs{i, col_new} = data(pos_indexes, :);  % Store parsed trials
        else
            pos_norm_PDpe_NTs{i, col_new} = [];  % Handle empty case
        end
        
        % Extract trials based on neg_norm_PDpe indexes and store in corresponding column
        if ~isempty(neg_indexes)
            neg_norm_PDpe_NTs{i, col_new} = data(neg_indexes, :);  
        else
            neg_norm_PDpe_NTs{i, col_new} = [];  % Handle empty case
        end
    end
end

%% concatenate across PD patients AUCs for each Neurotransmitter
for col = 1:3
    pos_norm_PDpe_NTs_concat{1, col} = vertcat(pos_norm_PDpe_NTs{:, col});
    neg_norm_PDpe_NTs_concat{1, col} = vertcat(neg_norm_PDpe_NTs{:, col});
end


for col = 1:3
    pos_norm_PDpe_NTs_mu{1, col} = mean(pos_norm_PDpe_NTs_concat{1, col}, 1, 'omitnan');
    pos_norm_PDpe_NTs_sem{1, col} = std(pos_norm_PDpe_NTs_concat{1, col}, 'omitnan') ./ sqrt(sum(~isnan(pos_norm_PDpe_NTs_concat{1, col}), 1));
    neg_norm_PDpe_NTs_mu{1, col} = mean(neg_norm_PDpe_NTs_concat{1, col}, 1, 'omitnan');
    neg_norm_PDpe_NTs_sem{1, col} = std(neg_norm_PDpe_NTs_concat{1, col}, 'omitnan') ./ sqrt(sum(~isnan(neg_norm_PDpe_NTs_concat{1, col}), 1));
    
end

%% ET
pos_norm_ETpe_NTs = cell(size(patients_behavior_ET, 1), 3);
neg_norm_ETpe_NTs = cell(size(patients_behavior_ET, 1), 3);
pos_norm_ETpe = pos_normpe(1:6,1);
neg_norm_ETpe = neg_normpe(1:6,1);
% Loop through each subject (row)
for i = 1:size(patients_behavior_ET, 1)
    % Get the index values for pos_norm_ETpe and neg_norm_ETpe for this subject
    pos_indexes = pos_norm_ETpe{i};
    neg_indexes = neg_norm_ETpe{i};
    
    for col = task_event
        data = patients_behavior_ET{i, col}(:, 1:71);
        

        col_new = col - 7;  
        
        if ~isempty(pos_indexes)
            pos_norm_ETpe_NTs{i, col_new} = data(pos_indexes, :); 
        else
            pos_norm_ETpe_NTs{i, col_new} = [];  
        end
        
        % Extract trials based on neg_norm_ETpe indexes and store in corresponding column
        if ~isempty(neg_indexes)
            neg_norm_ETpe_NTs{i, col_new} = data(neg_indexes, :);  
        else
            neg_norm_ETpe_NTs{i, col_new} = [];  
        end
    end
end
%% 
 for col = 1:3
    pos_norm_ETpe_NTs_concat{1, col} = vertcat(pos_norm_ETpe_NTs{:, col});
    neg_norm_ETpe_NTs_concat{1, col} = vertcat(neg_norm_ETpe_NTs{:, col});
 end

for col = 1:3
    pos_norm_ETpe_NTs_mu{1, col} = mean(pos_norm_ETpe_NTs_concat{1, col}, 1, 'omitnan');
    pos_norm_ETpe_NTs_sem{1, col} = std(pos_norm_ETpe_NTs_concat{1, col}, 'omitnan') ./ sqrt(sum(~isnan(pos_norm_ETpe_NTs_concat{1, col}), 1));
    neg_norm_ETpe_NTs_mu{1, col} = mean(neg_norm_ETpe_NTs_concat{1, col}, 1, 'omitnan');
    neg_norm_ETpe_NTs_sem{1, col} = std(neg_norm_ETpe_NTs_concat{1, col}, 'omitnan') ./ sqrt(sum(~isnan(neg_norm_ETpe_NTs_concat{1, col}), 1));
end


%% PD AUCs by trial
patients_behavior_PD_AUC(:,1:7) = patients_behavior_PD(:,1:7);

for subject = 1:length(patients_behavior_PD(:,1))
    for col= task_event
        for row= 1:length(patients_behavior_PD{subject,col}(:,1))
patients_behavior_PD_AUC{subject,col}(row,1) = trapz(patients_behavior_PD{subject, col}(row,analyze_window));
        end
    end
end
%% ET AUCs by trial
patients_behavior_ET_AUC(:,1:7) = patients_behavior_ET(:,1:7);

for subject = 1:length(patients_behavior_ET(:,1))
    for col= task_event
        for row= 1:length(patients_behavior_ET{subject,col}(:,1))
patients_behavior_ET_AUC{subject,col}(row,1) = trapz(patients_behavior_ET{subject, col}(row,analyze_window));
        end
    end
end

%% density curve Distribution for offers
[fH_ET, xiH_ET] = ksdensity(ET_caudate_high_offers);
[fL_ET, xiL_ET] = ksdensity(ET_caudate_low_offers);
[fH_PD, xiH_PD] = ksdensity(PD_caudate_high_offers);
[fL_PD, xiL_PD] = ksdensity(PD_caudate_low_offers);

%% Emotional Ratings ET

bar_start = -5;
% Initialize arrays to store the average and SEM for each emotional rating (1 to 9)
num_ratings = 9;
average_norm_pe = nan(1, num_ratings);
sem_norm_pe = nan(1, num_ratings);

% List of selected subjects
selected_subjects = find(strcmp(patients_table_behavior(:,2), "ET"))';

% Loop through each emotional rating (1 to 9)
for rating = 1:num_ratings
    % Collect all norm_pe values that correspond to the current emotional rating
    norm_pe_values_for_rating = [];

    % Loop through only the selected subjects
    for i = selected_subjects
        % Extract the combined norm_pe and emotional ratings for the subject
        combined_trials = norm_pe_with_emotion{i};

        matching_indices = combined_trials(:, 2) == rating;
        norm_pe_values_for_rating = [norm_pe_values_for_rating; combined_trials(matching_indices, 1)];
    end

    average_norm_pe_ET(rating) = mean(norm_pe_values_for_rating, 'omitnan');
    sem_norm_pe_ET(rating) = std(norm_pe_values_for_rating, 'omitnan') / sqrt(sum(~isnan(norm_pe_values_for_rating)));
end

%% Emotional Ratings PD

bar_start = -5;

% Initialize arrays to store the average and SEM for each emotional rating (1 to 9)
num_ratings = 9;
average_norm_pe = nan(1, num_ratings);
sem_norm_pe = nan(1, num_ratings);

% List of selected subjects
selected_subjects = find(strcmp(patients_table_behavior(:,2), "PD"))';

% Loop through each emotional rating (1 to 9)
for rating = 1:num_ratings
    % Collect all norm_pe values that correspond to the current emotional rating
    norm_pe_values_for_rating = [];

    % Loop through only the selected subjects
    for i = selected_subjects
        % Extract the combined norm_pe and emotional ratings for the subject
        combined_trials = norm_pe_with_emotion{i};

        % Get the norm_pe values where the emotional rating matches the current rating
        matching_indices = combined_trials(:, 2) == rating;
        norm_pe_values_for_rating = [norm_pe_values_for_rating; combined_trials(matching_indices, 1)];
    end

    % Compute the average and SEM for the norm_pe values corresponding to this rating
    average_norm_pe_PD(rating) = mean(norm_pe_values_for_rating, 'omitnan');
    sem_norm_pe_PD(rating) = std(norm_pe_values_for_rating, 'omitnan') / sqrt(sum(~isnan(norm_pe_values_for_rating)));
end

%% ET LME table

rowIndex = 1; % To keep track of the row index for patients_table, each row will be a trial for every patient.


% Loop through each subject
for i = 1:length(patients_behavior_ET_AUC(:,1))
    % Loop through each neurotransmitter (columns 8:10)
        % Loop through each trial for the current patient
        for trial = 1:length(pos_norm_ETpe{i, 1}) % Assuming pos_norm_ETpe is defined
            
            patients_table_ET_LME{rowIndex, 1} = i; %subject ID     
            patients_table_ET_LME{rowIndex, 2} = 1; % Patient type ET = 1
            patients_table_ET_LME{rowIndex, 3} = 1;  % pos PE = 1
            patients_table_ET_LME{rowIndex, 4} = patients_behavior_ET_AUC{i, 8}(pos_norm_ETpe{i, 1}(trial)); % Column 8 AUC value
            patients_table_ET_LME{rowIndex, 5} = patients_behavior_ET_AUC{i, 9}(pos_norm_ETpe{i, 1}(trial)); % Column 9 AUC value
            patients_table_ET_LME{rowIndex, 6} = patients_behavior_ET_AUC{i, 10}(pos_norm_ETpe{i, 1}(trial)); % Column 10 AUC value          
            patients_table_ET_LME{rowIndex, 7} = pos_norm_ETpe{i, 1}(trial); % trial index. out of 30 trials          
            patients_table_ET_LME{rowIndex, 8} = norm_ET_pe(i,pos_norm_ETpe{i, 1}(trial)); % PE value per trial          
            patients_table_ET_LME{rowIndex, 9} = (patients_behavior_ET_AUC{i, 7}(pos_norm_ETpe{i, 1}(trial),2) - patients_behavior_ET_AUC{i, 7}(pos_norm_ETpe{i, 1}(trial),1))-4; %%%reaction time from offer to submit
            patients_table_ET_LME{rowIndex, 10} = patients_behavior_ET_AUC{i, 6}(pos_norm_ETpe{i, 1}(trial),1); %%emotional rating
            patients_table_ET_LME{rowIndex, 11} = patients_behavior_ET_AUC{i, 5}(pos_norm_ETpe{i, 1}(trial),2); %%accept/reject
            
            
            % Increment the row index
            rowIndex = rowIndex + 1;
        end



        for trial = 1:length(neg_norm_ETpe{i, 1})
            
            patients_table_ET_LME{rowIndex, 1}  = i; %subject ID     
            patients_table_ET_LME{rowIndex, 2} = 1; % Patient type ET = 1
            patients_table_ET_LME{rowIndex, 3} = 2; % neg PE = 2
            patients_table_ET_LME{rowIndex, 4} = patients_behavior_ET_AUC{i, 8}(neg_norm_ETpe{i, 1}(trial)); % Column 8 AUC value
            patients_table_ET_LME{rowIndex, 5} = patients_behavior_ET_AUC{i, 9}(neg_norm_ETpe{i, 1}(trial)); % Column 9 AUC value
            patients_table_ET_LME{rowIndex, 6} = patients_behavior_ET_AUC{i, 10}(neg_norm_ETpe{i, 1}(trial)); % Column 10 AUC value
            patients_table_ET_LME{rowIndex, 7}  = neg_norm_ETpe{i, 1}(trial); % trial index. out of 30 trials          
            patients_table_ET_LME{rowIndex, 8}  = norm_ET_pe(i,neg_norm_ETpe{i, 1}(trial)); % PE value per trial          
            patients_table_ET_LME{rowIndex, 9}  = (patients_behavior_ET_AUC{i, 7}(neg_norm_ETpe{i, 1}(trial),2) - patients_behavior_ET_AUC{i, 7}(neg_norm_ETpe{i, 1}(trial),1))-4; %%%reaction time from offer to submit
            patients_table_ET_LME{rowIndex, 10} = patients_behavior_ET_AUC{i, 6}(neg_norm_ETpe{i, 1}(trial),1); %%emotional rating
            patients_table_ET_LME{rowIndex, 11} = patients_behavior_ET_AUC{i, 5}(neg_norm_ETpe{i, 1}(trial),2); %%accept/reject

            rowIndex = rowIndex + 1;

        end
end

% Convert the cell array to a table for better usability
patients_table_LME_ET_var = cell2table(patients_table_ET_LME, 'VariableNames', { 'subjectNum','PatientType', 'PEtype', 'AUCDA', 'AUCNE', 'AUC5HT',  'IdxTrial', 'PEvalue',...
    'reactionTime', 'emoRating',  'offerDecision', });

%% PD LME table

PD_subjectNum = i;
rowIndex = 1; % To keep track of the row index for patients_table, each row will be a trial for every patient.


% Loop through each subject
for i = 1:length(patients_behavior_PD_AUC(:,1))
    % Loop through each neurotransmitter (columns 8:10)
        % Loop through each trial for the current patient
        for trial = 1:length(pos_norm_PDpe{i, 1}) % Assuming pos_norm_PDpe is defined
            
            patients_table_PD_LME{rowIndex, 1} = PD_subjectNum + i; %subject ID     
            patients_table_PD_LME{rowIndex, 2} = 2; % Patient type PD = 2
            patients_table_PD_LME{rowIndex, 3} = 1;  % pos PE = 1
            patients_table_PD_LME{rowIndex, 4} = patients_behavior_PD_AUC{i, 8}(pos_norm_PDpe{i, 1}(trial)); % Column 8 AUC value
            patients_table_PD_LME{rowIndex, 5} = patients_behavior_PD_AUC{i, 9}(pos_norm_PDpe{i, 1}(trial)); % Column 9 AUC value
            patients_table_PD_LME{rowIndex, 6} = patients_behavior_PD_AUC{i, 10}(pos_norm_PDpe{i, 1}(trial)); % Column 10 AUC value          
            patients_table_PD_LME{rowIndex, 7} = pos_norm_PDpe{i, 1}(trial); % trial index. out of 30 trials          
            patients_table_PD_LME{rowIndex, 8} = norm_PD_pe(i,pos_norm_PDpe{i, 1}(trial)); % PE value per trial          
            patients_table_PD_LME{rowIndex, 9} = (patients_behavior_PD_AUC{i, 7}(pos_norm_PDpe{i, 1}(trial),2) - patients_behavior_PD_AUC{i, 7}(pos_norm_PDpe{i, 1}(trial),1))-4; %%%reaction time from offer to submit
            patients_table_PD_LME{rowIndex, 10} = patients_behavior_PD_AUC{i, 6}(pos_norm_PDpe{i, 1}(trial),1); %%emotional rating
            patients_table_PD_LME{rowIndex, 11} = patients_behavior_PD_AUC{i, 5}(pos_norm_PDpe{i, 1}(trial),2); %%accept/reject
            
            
            % Increment the row index
            rowIndex = rowIndex + 1;
        end



        for trial = 1:length(neg_norm_PDpe{i, 1})
            
            patients_table_PD_LME{rowIndex, 1}  = PD_subjectNum + i; %subject ID   
            patients_table_PD_LME{rowIndex, 2} = 2; % Patient type PD = 2
            patients_table_PD_LME{rowIndex, 3} = 2;                            % neg PE = 2
            patients_table_PD_LME{rowIndex, 4} = patients_behavior_PD_AUC{i, 8}(neg_norm_PDpe{i, 1}(trial)); % Column 8 AUC value
            patients_table_PD_LME{rowIndex, 5} = patients_behavior_PD_AUC{i, 9}(neg_norm_PDpe{i, 1}(trial)); % Column 9 AUC value
            patients_table_PD_LME{rowIndex, 6} = patients_behavior_PD_AUC{i, 10}(neg_norm_PDpe{i, 1}(trial)); % Column 10 AUC value
            patients_table_PD_LME{rowIndex, 7}  = neg_norm_PDpe{i, 1}(trial); % trial index. out of 30 trials          
            patients_table_PD_LME{rowIndex, 8}  = norm_PD_pe(i,neg_norm_PDpe{i, 1}(trial)); % PE value per trial          
            patients_table_PD_LME{rowIndex, 9}  = (patients_behavior_PD_AUC{i, 7}(neg_norm_PDpe{i, 1}(trial),2) - patients_behavior_PD_AUC{i, 7}(neg_norm_PDpe{i, 1}(trial),1))-4; %%%reaction time from offer to submit
            patients_table_PD_LME{rowIndex, 10} = patients_behavior_PD_AUC{i, 6}(neg_norm_PDpe{i, 1}(trial),1); %%emotional rating
            patients_table_PD_LME{rowIndex, 11} = patients_behavior_PD_AUC{i, 5}(neg_norm_PDpe{i, 1}(trial),2); %%accept/reject

            rowIndex = rowIndex + 1;

        end
end

% Convert the cell array to a table for better usability
patients_table_LME_PD_var = cell2table(patients_table_PD_LME, 'VariableNames', {'subjectNum','PatientType', 'PEtype', 'AUCDA', 'AUCNE', 'AUC5HT',  'IdxTrial', 'PEvalue',...
    'reactionTime', 'emoRating',  'offerDecision', });

%% combined ET and PD tables
patients_table_LME_ET_PD_combined = [patients_table_LME_ET_var; patients_table_LME_PD_var];

%% Save .csv for R studio statistics and additional figure panels 

writetable(patients_table_LME_ET_PD_combined, full_path);  % Save the table

%% organize rows by trial # per subject
nET     = 6;
nPD     = 13;
ET_data = cell2mat(patients_table_ET_LME);
PD_data = cell2mat(patients_table_PD_LME);

%%% 
for n = 1:nET
    dat_tmp = ET_data((((n-1)*30)+1):(n*30),:);
    dat_tmp = sortrows(dat_tmp,7);
    ET_data((((n-1)*30)+1):(n*30),:) = dat_tmp;
end
for n = 1:nPD
    dat_tmp = PD_data((((n-1)*30)+1):(n*30),:);
    dat_tmp = sortrows(dat_tmp,7);
    PD_data((((n-1)*30)+1):(n*30),:) = dat_tmp;
end

UG_data = [ET_data ; PD_data];
%% find mean AUC for DA, NA and 5HT for Pos/Neg NPE per subject
nET     = 6;
nPD     = 13;

ET_data = UG_data(UG_data(:,2) == 1, :);

ET_data_PE = zeros(nET,6);
for n = 1:nET
    ET_data_PE(n,1) = mean(ET_data(ET_data((((n-1)*30)+1):(n*30),3)==1,4),'omitnan');
    ET_data_PE(n,2) = mean(ET_data(ET_data((((n-1)*30)+1):(n*30),3)==2,4),'omitnan');
    ET_data_PE(n,3) = mean(ET_data(ET_data((((n-1)*30)+1):(n*30),3)==1,5),'omitnan');
    ET_data_PE(n,4) = mean(ET_data(ET_data((((n-1)*30)+1):(n*30),3)==2,5),'omitnan');
    ET_data_PE(n,5) = mean(ET_data(ET_data((((n-1)*30)+1):(n*30),3)==1,6),'omitnan');
    ET_data_PE(n,6) = mean(ET_data(ET_data((((n-1)*30)+1):(n*30),3)==2,6),'omitnan');

end

PD_data = UG_data(UG_data(:,2) == 2, :);

PD_data_PE = zeros(nPD,6);
for n = 1:nPD
    PD_data_PE(n,1) = mean(PD_data(PD_data((((n-1)*30)+1):(n*30),3)==1,4),'omitnan');
    PD_data_PE(n,2) = mean(PD_data(PD_data((((n-1)*30)+1):(n*30),3)==2,4),'omitnan');
    PD_data_PE(n,3) = mean(PD_data(PD_data((((n-1)*30)+1):(n*30),3)==1,5),'omitnan');
    PD_data_PE(n,4) = mean(PD_data(PD_data((((n-1)*30)+1):(n*30),3)==2,5),'omitnan');
    PD_data_PE(n,5) = mean(PD_data(PD_data((((n-1)*30)+1):(n*30),3)==1,6),'omitnan');
    PD_data_PE(n,6) = mean(PD_data(PD_data((((n-1)*30)+1):(n*30),3)==2,6),'omitnan');

end

%% Compute SVD based on mean AUC data

PE_svd_dat = [ET_data_PE ; PD_data_PE];
[u,s,v]    = svd(PE_svd_dat,0);
vt         = v';

%% Variance Explained by Each SV Supplementary Table 6
singular_vals = diag(s);                  
explained_var = singular_vals.^2;        
explained_var_ratio = explained_var / sum(explained_var);  
explained_var_percent = 100 * explained_var_ratio;

%% Logistic classifier analysis to predict ET/PD class from SVs
niters = 10000;

model_idx_single(1:6,1) = 1:6; 
model_idx_single(1:6,2) = zeros(6,1);

model_idx_paired(1:5,1)   = ones(5,1);   model_idx_paired(1:5,2)   = 2:6;
model_idx_paired(6:9,1)   = 2*ones(4,1); model_idx_paired(6:9,2)   = 3:6;
model_idx_paired(10:12,1) = 3*ones(3,1); model_idx_paired(10:12,2) = 4:6;
model_idx_paired(13:14,1) = 4*ones(2,1); model_idx_paired(13:14,2) = 5:6;
model_idx_paired(15,1)    = 5;           model_idx_paired(15,2)    = 6;

model_idx = [model_idx_single ; model_idx_paired];
AUC_perm  = zeros(size(model_idx,1),niters); 

svd_reg_results = zeros(size(model_idx,1),2);

% do SVD
PE_svd_dat     = [ET_data_PE ; PD_data_PE];
y              = [ones(nET,1) ; zeros(nPD,1)];  % ET = 1 , PD = 0;
idx_tmp        = 1:nParticipants;
PE_svd_dat_tmp = PE_svd_dat(idx_tmp,:);
[u,s,v]        = svd(PE_svd_dat_tmp,0);

% do svd regression and permutation tests for each model
for m = 1:size(model_idx,1)

    idxs_tmp = model_idx(m,:);
    idxs     = idxs_tmp(idxs_tmp ~= 0);
    x        = u(:,idxs);

    disp(['Permutation test model ' num2str(m) ' / ' num2str(size(model_idx,1)) ' , SV(s): ' num2str(idxs) ])

    Log_ET_PD        = fitglm(x,y,'Distribution','binomial','Link','logit','LikelihoodPenalty','jeffreys-prior');
    score_log        = Log_ET_PD.Fitted.Probability;
    [~,~,~,AUC(m,1)] = perfcurve(y,score_log,1);

    for i = 1:niters    
        if mod(i,10000) == 0
            disp(['iteration ' num2str(i) ' / ' num2str(niters)])
        end
        perm_idx              = randperm(numel(y));
        ytmp                  = y(perm_idx);
        Log_ET_PD_iter        = fitglm(x,ytmp,'Distribution','binomial','Link','logit','LikelihoodPenalty','jeffreys-prior');
        score_log             = Log_ET_PD_iter.Fitted.Probability;
        [~,~,~,AUC_perm(m,i)] = perfcurve(ytmp,score_log,1);
    end

    svd_reg_results(m,1) = AUC(m,1);
    svd_reg_results(m,2) = sum(AUC_perm(m,:)>=AUC(m,1))/niters;
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Plotting
ET_color = [80, 233, 145]/ 255;
PD_color = [155, 25, 245]/ 255;

mkralpha = 0.7;
lineW = 3;
%%  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Plotting Figure 1b Low Offer Distribution

dis_fig = figure;
hold on
plot(xiL_ET, fL_ET, 'LineWidth', 2, 'Color', ET_color);  
plot(xiL_PD, fL_PD, 'LineWidth', 2, 'Color', PD_color);  
xlabel('Offers ($)');
ylabel('Sampling Frequency');
et_dis_title = 'Distribution Curve Low Offers';
title(et_dis_title);
    ax = gca;
        ax.FontSize = 20;
xlim([0,16])

hold off
%%  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Plotting Figure 1c High Offer Distribution
dis_fig = figure;
hold on
plot(xiH_ET, fH_ET, 'LineWidth', 2, 'Color', ET_color);  
plot(xiH_PD, fH_PD, 'LineWidth', 2, 'Color', PD_color);  
xlabel('Offers ($)');
ylabel('Sampling Frequency');
et_dis_title = 'Distribution Curve High Offers';
title(et_dis_title);
    ax = gca;
        ax.FontSize = 20;
xlim([0,16])

hold off
%%  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Plotting Figure 1e NPE Distribution


all_values_ET = norm_ET_pe(:);
all_values_PD = norm_PD_pe(:);

[f_ET, xi_ET] = ksdensity(all_values_ET); % KDE for ET
[f_PD, xi_PD] = ksdensity(all_values_PD); % KDE for PD

fig = figure; hold on;
plot(xi_ET, f_ET, 'Color', ET_color, 'LineWidth', 2); % ET curve (blue)
plot(xi_PD, f_PD, 'Color', PD_color ,'LineWidth', 2); % PD curve (red)
hold off;
yticks(0:0.05:0.1)
xlabel('NPE Value');
ylabel('Sampling Frequency');
title('Distribution of NPEs');
legend('ET', 'PD');
    ax = gca;
    ax.FontSize = 20
grid on;




%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Plotting  Figure 2a ET Positive NPE
ylimit_full_trace = [-.49, .46];
xlimit  = [3,7];
fontS = 25;
fs = 10;
time = (1:length(pos_norm_ETpe_NTs_mu{1,1}))/fs;
fig = figure;
hold on
    plot(time,pos_norm_ETpe_NTs_mu{1,1},'LineWidth',lineW,'color','k');    
    plot(time,pos_norm_ETpe_NTs_mu{1,2},'LineWidth',lineW,'color',"c"); 
    plot(time,pos_norm_ETpe_NTs_mu{1,3},'LineWidth',lineW,'color',"m"); 

    shadedErrorBar(time,pos_norm_ETpe_NTs_mu{1,1},pos_norm_ETpe_NTs_sem{1,1},'lineProps','k');
    shadedErrorBar(time,pos_norm_ETpe_NTs_mu{1,2},pos_norm_ETpe_NTs_sem{1,2},'lineProps', 'c');
    shadedErrorBar(time,pos_norm_ETpe_NTs_mu{1,3},pos_norm_ETpe_NTs_sem{1,3},'lineProps', 'm');

ylabel('NT [z]');
    xlabel('Time (s)');
ttl_onset = 3;

ylim(ylimit_full_trace)
xticks(0:1:7)
xticklabels({'-3', '-2', '-1', '0', '1', '2', '3', '4'})
xlim(xlimit)
    ax = gca;
    ax.FontSize = fontS;
 
hold off

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Plotting  Figure 2b ET Negative NPE

fig = figure;
hold on
    plot(time,neg_norm_ETpe_NTs_mu{1,1},'LineWidth',lineW,'color','k');    
    plot(time,neg_norm_ETpe_NTs_mu{1,2},'LineWidth',lineW,'color',"c"); 
    plot(time,neg_norm_ETpe_NTs_mu{1,3},'LineWidth',lineW,'color',"m"); 

    shadedErrorBar(time,neg_norm_ETpe_NTs_mu{1,1},neg_norm_ETpe_NTs_sem{1,1},'lineProps','k');
    shadedErrorBar(time,neg_norm_ETpe_NTs_mu{1,2},neg_norm_ETpe_NTs_sem{1,2},'lineProps', 'c');
    shadedErrorBar(time,neg_norm_ETpe_NTs_mu{1,3},neg_norm_ETpe_NTs_sem{1,3},'lineProps', 'm');

ylabel('NT [z]');
    xlabel('Time (s)');
ttl_onset = 3;
ylim(ylimit_full_trace)
xticks(0:1:7)
xticklabels({'-3', '-2', '-1', '0', '1', '2', '3', '4'})
xlim(xlimit)
    ax = gca;
    ax.FontSize = fontS;

hold off

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Plotting  Figure 2d PD Positive NPE
% postive PEs PD mean full snippet

time = (1:length(pos_norm_PDpe_NTs_mu{1,1}))/fs;
fig = figure;
hold on
    plot(time,pos_norm_PDpe_NTs_mu{1,1},'LineWidth',lineW,'color','k');    
    plot(time,pos_norm_PDpe_NTs_mu{1,2},'LineWidth',lineW,'color',"c"); 
    plot(time,pos_norm_PDpe_NTs_mu{1,3},'LineWidth',lineW,'color',"m"); 

    shadedErrorBar(time,pos_norm_PDpe_NTs_mu{1,1},pos_norm_PDpe_NTs_sem{1,1},'lineProps','k');
    shadedErrorBar(time,pos_norm_PDpe_NTs_mu{1,2},pos_norm_PDpe_NTs_sem{1,2},'lineProps', 'c');
    shadedErrorBar(time,pos_norm_PDpe_NTs_mu{1,3},pos_norm_PDpe_NTs_sem{1,3},'lineProps', 'm');

ylabel('NT [z]');
    xlabel('Time (s)');
ttl_onsPD = 3;
ylim(ylimit_full_trace)
xticks(0:1:7)
xticklabels({'-3', '-2', '-1', '0', '1', '2', '3', '4'})
xlim(xlimit)
    ax = gca;
    ax.FontSize = fontS;

hold off

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Plotting  Figure 2e PD Negative NPE

fig = figure;
hold on
    plot(time,neg_norm_PDpe_NTs_mu{1,1},'LineWidth',lineW,'color','k');    
    plot(time,neg_norm_PDpe_NTs_mu{1,2},'LineWidth',lineW,'color',"c"); 
    plot(time,neg_norm_PDpe_NTs_mu{1,3},'LineWidth',lineW,'color',"m"); 

    shadedErrorBar(time,neg_norm_PDpe_NTs_mu{1,1},neg_norm_PDpe_NTs_sem{1,1},'lineProps','k');
    shadedErrorBar(time,neg_norm_PDpe_NTs_mu{1,2},neg_norm_PDpe_NTs_sem{1,2},'lineProps', 'c');
    shadedErrorBar(time,neg_norm_PDpe_NTs_mu{1,3},neg_norm_PDpe_NTs_sem{1,3},'lineProps', 'm');

ylabel('NT [z]');
    xlabel('Time (s)');
ttl_onsPD = 3;

ylim(ylimit_full_trace)

xticks(0:1:7)
xticklabels({'-3', '-2', '-1', '0', '1', '2', '3', '4'})
xlim(xlimit)

    ax = gca;
    ax.FontSize = fontS;


hold off

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Plotting  Figure 2g SVD Score

fig = figure; 
hold on
for i = 1:nParticipants
    if i == 1 || i == 4 || i == 5 || i == 6
        bar([zeros(1,i-1) u(i,5) zeros(1,nParticipants-i)],'FaceColor',ET_color,'EdgeColor','k','FaceAlpha',mkralpha)
    elseif i == 2 || i == 3
        bar([zeros(1,i-1) u(i,5) zeros(1,nParticipants-i)],'FaceColor',ET_color,'EdgeColor','k','FaceAlpha',mkralpha)
    elseif i == 16 || i == 18 || i == 19
        bar([zeros(1,i-1) u(i,5) zeros(1,nParticipants-i)],'FaceColor',PD_color,'EdgeColor','k', 'FaceAlpha',mkralpha)
    else
        bar([zeros(1,i-1) u(i,5) zeros(1,nParticipants-i)],'FaceColor',PD_color,'EdgeColor','k', 'FaceAlpha',mkralpha)
    end
end

set(gca,'FontSize',24)
ylabel('SV 5 Score')
xticklabels([])
xlabel('Patient')

hold off

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Plotting  Figure 2h SVD Weight

fig = figure; 
hold on
bar([vt(5,1:2) , zeros(1,4)],'FaceColor','black')
bar([zeros(1,2) , vt(5,3:4) , zeros(1,2)],'FaceColor','cyan')
bar([zeros(1,4) , vt(5,5:6)],'FaceColor','magenta')
xticks(1:6)
xticklabels({'Pos','Neg','Pos','Neg','Pos','Neg'})
ylabel('SV 5 weight')
set(gca,'FontSize',24)

hold off

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Plotting  Figure 2i SVD 3D hyperplane defined by right singular vectors 1, 2, and 5

mkrsize = 30;
u_idxs = [1 2 5];
y = [ones(nET,1) ; -1*ones(nPD,1)];     
x = [u(:,u_idxs) , ones(size(u,1),1)];
w = (x' * x) \ (x' * y);
w1 = w(1); w2 = w(2); w3 = w(3); b = w(4);
[x1grid,x2grid] = meshgrid(linspace(min(x(:,1)),max(x(:,1)),30),linspace(min(x(:,2)),max(x(:,2)),30));
x3grid = -((w1 * x1grid) + (w2 * x2grid) + b) / w3;


fig = figure
set(gcf,'color','w')
hold on
for i = 1:6
    hold on
        scatter3(u(i,1),u(i,2),u(i,5),mkrsize,'filled','MarkerFaceColor',ET_color,'MarkerEdgeColor','k','MarkerFaceAlpha',mkralpha,'HandleVisibility','off')
end

for i = 1:13
    hold on
        scatter3(u(i+nET,1),u(i+nET,2),u(i+nET,5),mkrsize,'filled','MarkerFaceColor',PD_color,'MarkerEdgeColor','k','MarkerFaceAlpha',mkralpha,'HandleVisibility','off')

end
C = zeros(size(x1grid,1),size(x1grid,1),3);
C(:,:,1) = 0.43*ones(size(x1grid,1)); C(:,:,2) = 0.42*ones(size(x1grid,1)) ; C(:,:,3) = 0.64*ones(size(x1grid,1));
surf(x1grid,x2grid,x3grid,C,'FaceAlpha',0.5,'EdgeColor','none')
xlabel('SV1')
ylabel('SV2')
zlabel('SV5')
fig.Units = 'inches';
fig.Position(3) = 2.3;  
fig.Position(4) = 1.8562; 
ax = gca;
ax.FontName = 'Arial';
ax.FontSize = 8; 
ax.LineWidth = 1; 
ax.Title.FontName = 'Arial';
ax.Title.FontSize = 9; 
view(38.3937,11.5115)
grid on

hold off

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Plotting  Figure 2j Mean DA Pos/Neg per subject


set(gcf,'color','w')
fs = 20;
ylims = [-6 6];
xlims = [-12 12];
mkrsize = 300;
mkralpha = 0.7;

fig = figure;
hold on
line([xlims(1) xlims(2)],[0 0],'LineWidth',1,'Color',[.5 .5 .5],'HandleVisibility','off')
line([0 0],[ylims(1) ylims(2)],'LineWidth',1,'Color',[.5 .5 .5],'HandleVisibility','off')
for i = 1:nET
    hold on
        scatter(ET_data_PE(i,1),ET_data_PE(i,2),mkrsize,'filled','MarkerFaceColor',ET_color,'MarkerEdgeColor','k','MarkerFaceAlpha',mkralpha,'HandleVisibility','off')
    end
for i = 1:nPD
    hold on
        scatter(PD_data_PE(i,1),PD_data_PE(i,2),mkrsize,'filled','MarkerFaceColor',PD_color,'MarkerEdgeColor','k','MarkerFaceAlpha',mkralpha,'HandleVisibility','off')
    end

set(gca,'FontSize',fs)
xlabel('Mean DA AUC, Positive NPE')
ylabel('Mean DA AUC, Negative NPE')
axis square
ylim(ylims)
xlim(xlims)

hold off
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Plotting  Figure 2k Mean NA Pos/Neg per subject

fig = figure;
hold on
line([xlims(1) xlims(2)],[0 0],'LineWidth',1,'Color',[.5 .5 .5],'HandleVisibility','off')
line([0 0],[ylims(1) ylims(2)],'LineWidth',1,'Color',[.5 .5 .5],'HandleVisibility','off')
for i = 1:6
    hold on
        scatter(ET_data_PE(i,3),ET_data_PE(i,4),mkrsize,'filled','MarkerFaceColor',ET_color,'MarkerEdgeColor','k','MarkerFaceAlpha',mkralpha,'HandleVisibility','off')
end
for i = 1:13
    hold on
        scatter(PD_data_PE(i,3),PD_data_PE(i,4),mkrsize,'filled','MarkerFaceColor',PD_color,'MarkerEdgeColor','k','MarkerFaceAlpha',mkralpha,'HandleVisibility','off')

end
set(gca,'FontSize',fs)
xlabel('Mean NA AUC, Positive NPE')
ylabel('Mean NA AUC, Negative NPE')
axis square
ylim(ylims)
xlim(xlims)

hold off
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Plotting  Figure 2l Mean 5HT Pos/Neg per subject

fig = figure; 
hold on
line([xlims(1) xlims(2)],[0 0],'LineWidth',1,'Color',[.5 .5 .5],'HandleVisibility','off')
line([0 0],[ylims(1) ylims(2)],'LineWidth',1,'Color',[.5 .5 .5],'HandleVisibility','off')
for i = 1:6
    hold on
        scatter(ET_data_PE(i,5),ET_data_PE(i,6),mkrsize,'filled','MarkerFaceColor',ET_color,'MarkerEdgeColor','k','MarkerFaceAlpha',mkralpha,'HandleVisibility','off')
end

for i = 1:13
    hold on
        scatter(PD_data_PE(i,5),PD_data_PE(i,6),mkrsize,'filled','MarkerFaceColor',PD_color,'MarkerEdgeColor','k','MarkerFaceAlpha',mkralpha,'HandleVisibility','off')

end
set(gca,'FontSize',fs)
xlabel('Mean 5HT AUC, Positive NPE')
ylabel('Mean 5HT AUC, Negative NPE')
axis square
ylim(ylims)
xlim(xlims)

hold off
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%Plotting Supplementary Figure 2a Norm Prediction Error

ax = gca;
ax.SortMethod = 'childorder';
fontsize_tmp = 20;
fig = figure;
hold on
for i = 1:nParticipants
    if i == 1 || i == 4 || i == 5 || i == 6
       plot(norm_pe(i,:),':','Color',ET_color,'LineWidth',2)
    elseif i == 2 || i == 3
       plot(norm_pe(i,:),'Color',ET_color,'LineWidth',2)
    elseif i == 16 || i == 18 || i == 19
       h= plot(norm_pe(i,:),':','Color',PD_color,'LineWidth',2)
                    uistack(h, 'bottom')

    else
       h= plot(norm_pe(i,:),'Color',PD_color,'LineWidth',2)
                    uistack(h, 'bottom')

    end
end
xline(15, '--', 'LineWidth', 2)
ylabel('Norm Prediction Error')
xlabel('Trial')
set(gca,'fontsize',fontsize_tmp)
hold off



%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%Plotting Supplementary Figure 2b Bayesian Observer Norm
fig = figure;
set(gcf,'color','w')
hold on
for i = 1:nParticipants
    if i == 1 || i == 4 || i == 5 || i == 6
       plot(norm(i,:), ':', 'Color', ET_color, 'LineWidth', 2);
    elseif i == 2 || i == 3
       plot(norm(i,:),'Color',ET_color,'LineWidth',2)
    elseif i == 16 || i == 18 || i == 19
       h= plot(norm(i,:),':','Color',PD_color,'LineWidth',2)
             uistack(h, 'bottom')
    else
       h= plot(norm(i,:),'Color',PD_color,'LineWidth',2)
                    uistack(h, 'bottom')

    end
end

xline(15, '--', 'LineWidth', 2)
ylim([1 16])
ylabel('Bayesian Observer Norm')
xlabel('Trial')
set(gca,'fontsize',fontsize_tmp)
hold off



%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Plotting Supplementary Figure 2e PD, NPEs per Emotional Rating
fig = figure;
hold on;

for i = 1:num_ratings
    bar_width = 0.6;
    x_bar = [i - bar_width/2, i - bar_width/2, i + bar_width/2, i + bar_width/2];
    y_bar = [bar_start, average_norm_pe_PD(i), average_norm_pe_PD(i), bar_start];

    fill(x_bar, y_bar, PD_color, 'EdgeColor', 'k','FaceAlpha',mkralpha); 
end

for i = 1:num_ratings
    errorbar(i, average_norm_pe_PD(i), sem_norm_pe_PD(i), 'k', 'LineWidth', 1.5, 'CapSize', 10);
end

xlabel('Emotional Rating');
ylabel('Norm Prediction Errors');
title('PD')
xticks(1:num_ratings);
xticklabels({'1', '2', '3', '4', '5', '6', '7', '8', '9'});
ylim([bar_start, 4.5]);
hold off;
set(gca,'fontsize',fontsize_tmp)



%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Plotting Supplementary Figure 2f ET, NPEs per Emotional Rating

fig = figure;
hold on;

for i = 1:num_ratings
    bar_width = 0.6;
    x_bar = [i - bar_width/2, i - bar_width/2, i + bar_width/2, i + bar_width/2];
    y_bar = [bar_start, average_norm_pe_ET(i), average_norm_pe_ET(i), bar_start];

    fill(x_bar, y_bar, ET_color, 'EdgeColor', 'k','FaceAlpha',mkralpha); 
end

for i = 1:num_ratings
    errorbar(i, average_norm_pe_ET(i), sem_norm_pe_ET(i), 'k', 'LineWidth', 1.5, 'CapSize', 10);
end

xlabel('Emotional Rating');
ylabel('Norm Prediction Errors');
title('ET')
xticks(1:num_ratings);
xticklabels({'1', '2', '3', '4', '5', '6', '7', '8', '9'});
ylim([bar_start, 4.5]);
hold off;
set(gca,'fontsize',fontsize_tmp)







