function autism_paper_nCC_nPL_MD_parrallel_Github(x)

% -----------------------------------------------------------------------
% This script was produced and tested by Chirag Mehra, for the work found in the manuscript: 
% Mehra et al., (2026): "EEG functional connectivity as a prognostic biomarker of adaptive function in autistic people"
% Please cite the most up to date version of the manuscript when using this script

% This script was designed for parallel processing of one participant's EEG
% data at a time. It can be edited to processs participant data in series.
% It calls upon scripts from the Brain Connectivity Toolbox.

% Input: 1-4 functional connectivity matricies per participant
% Output: row containing mean degree, normalised weighted clustering
% coefficient, normalised weighted pathlength for each of 1-4 EEG epochs.
% -----------------------------------------------------------------------

%set up paths
addpath('')
addpath('/.../Brain_connectivity_toolbox/')

MainDirectory = '...'; % define a main directory
netfiles = dir(fullfile(MainDirectory,'nets_ROIdata*'));

path2output = '...';

output = nan*ones(1,13);

%run loop

for i = 1:length(netfiles)

    %%
    disp(i)

    cd(netfiles(i).folder);
    data = load(netfiles(i).name);
    fnames = fieldnames(data);
    networks = data.(fnames{1});
    clear data
    clear fnames

    n_networks = size(networks,3); %number of networks in each file

    %make temporary output variables
    md_temp = nan*ones(1,4);
    cc_temp = nan*ones(1,4);
    pl_temp = nan*ones(1,4);

    for j = 1:n_networks %run each network for a participant one at a time

        %loop through networks per person
        connectivity_matrix = squeeze(networks(:,:,j)); %squeeze makes 3D into 2D by removing dimensions of length 1

        %clustering coefficient code
        cc_vector = clustering_coef_wu(connectivity_matrix);
        mean_cc = mean(cc_vector);

        %path length code
        temp_matrix = weight_conversion(connectivity_matrix, 'lengths'); %weights into distances by inversion method
        [D,~] = distance_wei(temp_matrix); %calculates the distance between all pairs of nodes
        [lambda,~,~,~,~] = charpath(D); %calculates the average distance between all pairs of nodes

        % calculate CC and PL for surrogate networks
        surrogate_matrix = nan*ones(500,2);

        for k = 1:500

            %make null models
            [W0,~] = null_model_und_sign(connectivity_matrix,5,0.5); %Random graphs with preserved weight, degree and
            %strength distributions

            % calculate CC for surrogate network
            CC_sug = clustering_coef_wu(W0);
            % mean CC across all nodes
            meanCC_sug = mean(CC_sug);
            %save CC into matrix
            surrogate_matrix(k,1) = meanCC_sug;

            %path length stuff
            temp_matrix_sug = weight_conversion(W0, 'lengths'); %weights into distances by inversion method
            [D_sug,~] = distance_wei(temp_matrix_sug); %calculates the distance between all pairs of nodes
            [lambda_sug,~,~,~,~] = charpath(D_sug,0,0); %calculates the average distance between all pairs of nodes

            % save surrogate PL into surrogate matrix
            surrogate_matrix(k,2) = lambda_sug;

        end

        nCC = (mean_cc/nanmean(surrogate_matrix(:,1)));
        nPL = (lambda/nanmean(surrogate_matrix(:,2)));

        cc_temp(1,j) = nCC;
        pl_temp(1,j) = nPL;

        %mean degree code
        [deg] = degrees_und(connectivity_matrix); %undirected
        md = mean(deg);

        md_temp(1,j) = md;

    end

    %create output variable
    pre_ID = erase(string(netfiles(i).name), ".mat");
    ID = str2double(erase(string(pre_ID), "nets_ROIdata "));
    output = [ID md_temp cc_temp pl_temp];

    %%
    %display data as a table
    Summary_data = array2table(output);

    %Name variables - THERE ARE 25 OF THESE
    Summary_data.Properties.VariableNames = {'Subject' ...
        'mean_degree_1' 'mean_degree_2' 'mean_degree_3' 'mean_degree_4'...
        'cc_1' 'cc_2' 'cc_3' 'cc_4'...
        'pl_1' 'pl_2' 'pl_3' 'pl_4'...
        };

    cd(path2output)
    save_name = join(['graph_measures_',num2str(ID)]);
    save(save_name, 'Summary_data')

end

exit