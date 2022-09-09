%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                              %
%                     DETECTION OF SEGMENTATION BOUNDARIES                     %
%                      AND PRODUCTION OF A MOTION PICTURE                      %  
%                  OF MULTI-MARKER OPTICAL MOTION CAPTURE DATA                 %
%                                                                              %
%                                    May 2019                                  %
%                                                                              %
%                          Juan Ignacio Mendoza Garay                          %
%                               doctoral student                               %
%                 Department of Music, Art and Culture Studies                 %
%                            University of Jyv?skyl?                           %
%                                                                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% This program has been tested with:
%   Matlab R2015a
%   Mocap Toolbox v1.5 

% The Mocap Toolbox can be downloaded from this webpage:
% https://www.jyu.fi/hytk/fi/laitokset/mutku/en/research/materials/mocaptoolbox

% However, the functions 'mcanimate' and 'mcplotframe' of the Mocap Toolbox 
% needed to be modified.
% The modified functions are provided with this script.

% For distances other than 'euclidean' the program requires the 'pdsit' function,
% included in Mathworks' Statistics and Machine Learning Toolbox.

% ==============================================================================
% Description:

% This program first computes segmentation boundaries in marker-based 
% motion-captured data.
% Then, segments can be rendered to stills or motion pictures separately.
% The input data files contain 3D marker trajectories in global coordinates and 
% can have extension 'tsv' or 'c3d'.

% The required files of motion-capture recordings are:
%   1) man_heaven_six.tsv
%   2) 14_20.c3d
%   3) 18_15.c3d
%   4) man_woman_lindy_hop.tsv
%   5) woman_disco.tsv

% (1), (4) and (5) are provided with this script.
% (2) and (3) can be downloaded from this webpage: http://mocap.cs.cmu.edu/search.php

% Two off-line methods and one on-line method for detection of segmentation 
% boundaries can be used, which have been developed using the principles 
% described by Foote (2000).

% ==============================================================================
% Instructions:

% The unmodified program will produce the results reported in the scholarly article,
% except computation time which depends on the computing machine used.
% To obtain different results, edit parameters marked with an arrow like this: <===
% To make further adjustments, edit parameters marked with an arrow like this: <---
% Note that each group of variables below "INFORMATION ABOUT DATA AND VARIABLES" 
% are for one motion-capture recording. Therefore, a new recording can be incorporated following the
% same logic.
 
% ==============================================================================
% Initialisation:

clc
clear
close all
restoredefaultpath

% ==============================================================================
% Declare paths:

  toolbox_path = ' '; % <--- folder where the Mocap toolbox is
     data_path = ' '; % <--- folder where data files are
  results_path = ' '; % <--- folder to save results, including still and motion pictures
     
addpath(genpath(toolbox_path))

%% -----------------------------------------------------------------------------

% DESCRIPTION OF RECORDINGS:
%   1 = JYU, a man performing alternated weaving of sticks and rest 
%   2 = CMU, subject 14, trial 20; jumping, twisting, raising knees, reaching up and down, reaching right and left, arm circles
%   3 = CMU, subjects 18 and 19, trial 15; Chicken Dance 2 times continuously (although only data of subject A/18 is used)
%   4 = JYU, a man and a woman dancing Lindy Hop style to the first 126.1 seconds of Gordon Webster's "I like you best of all", at approx. 167 bpm, total 350 beats (although only data of woman is used)
%   5 = JYU, a woman moving quasi-spontaneously to the first 108.5 seconds of The Bee Gees' "Stayin' Alive", at approx. 103.7 bpm, total 187 beats

         all_datainfo = 1:5;         % <=== select one or more motion-capture recordings
     all_resample_fps = 30;          % <=== frames per second for resampling (scalar or vector)
all_derivative_orders = 1;           % <=== select one or two derivative orders (0 = position, 1 = velocity, [0,1] = both)
      all_granularity = 'all';       % <=== select one or more granularity levels (1 = most coarse, greater number = finer granularity, or 'all')
  all_distance_metric = 'euclidean'; % <=== distance metric ('euclidean', 'cityblock', etc., see 'help pdist', put into cell for more than one: {'cityblock','euclidean'} )
           all_method = 0;           % <=== use on-line or off-line algorithm ( 0 = on-line, 1 = off-line as online, 2 = off-line alternative, [0,1,2] = all ; offline-as-online will give same result as online for some metrics i.e., euclidean)
       plotprocess_sw = 0;           % <=== 1 = plot online real-time process or offline distance matrix, 0 = don't plot them
           segplot_sw = 0;           % <=== plot each segmentation sequence in a separate figure (0 = no, 1 = yes)
            segnum_sw = 1;           % <=== display segment numbers in single-sequence segmentation plot (0 = no, 1 = yes)
plotsaveallresults_sw = 1;           % <=== save all results and a figure with all segmentation plots (0 = no, 1 = yes)
  segplot_save_format = 'tiff';      % <=== file format for saving segmentation plot figures ('fig' or 'tiff', {'fig','tiff'} for both, [] = don't save)   
      pic_save_format = [];          % <=== file format for saving pictures of markers ('fig' or 'tiff' for still, 'mpeg4' for video, {'fig','mpeg4'} or {'tiff','mpeg4'} for both or  {'fig', 'tiff','mpeg4'} for all, [] = don't make or save)) 
all_selected_segments = 'all';       % <=== select segments for pictures (e.g., [3:5] = segments 3 to 5, 'all' = all segments, [] = make pictures of un-segmented data without doing segmentation)
           
           economy_sw = 1;           % <--- 1 = don't compute segmentation processes if results are in memory, 0 = compute anyway, [] = don't compute segmentation
   video_trace_length = 0.2;         % <--- maximum trace length for moving picture (in seconds, 'full' = trace whole segment)
   preview_frame_time = [];          % <--- frames to preview (in seconds, [] = no preview)
        colour_switch = 0;           % <--- for segmentation plots, 0 = gray, 1 = colours
       select_monitor = 1;           % <--- for segmentation plots, monitor to display in
          dslabels_sw = 1;           % <--- for segmentation plots, 0 = use data file names for titles; 1 = use recording index
                                            
% ------------------------------------------------------------------------------

close all
clc
program_start_time = tic;
if iscell(all_distance_metric) == 0
    all_distance_metric = {all_distance_metric};
end

% initialise results record:
clear results
results.recording_filename = {0};
results.distance_measure = [];
results.numvars_headers = [{'recording'},{'n_markers'},{'fps'},{'derivative_order'},...
    {'novelty_width_f'},{'filter_kernel_width_f'},{'peaks_threshold'},{'novelty_width_f/filter_width_f'}];
results.numvars_values = nan(1,8);
results.nov_filt = {0};   % filtered novelty score
results.seg_bounds = {0}; % selected novelty peaks indexes (in frames)
results.times_headers = [{'data_time'},{'computation_time'},{'data_time/computation_time'}];
results.times_values = nan(1,3); % times (in seconds)
results_count = 0;

% parameters for plots:
oneplot_title_fontsize      = 22;
oneplot_ticklabels_fontsize = 18;
oneplot_axeslabels_fontsize = 22;
if colour_switch
    seg_boundaries_colour = [0.2 0.7 0.2];
    seg_curve_colour = [0 0 0.8];
    seg_labels_colour = seg_boundaries_colour - 0.2;
    seg_threshold_colour = seg_boundaries_colour + 0.1;
else
    seg_boundaries_colour = [1 1 1] * 0.6;
    seg_curve_colour = [0 0 0];
    cmap = gray(256);
    seg_labels_colour = seg_boundaries_colour - 0.5;
    seg_threshold_colour = seg_labels_colour;
end
monpos = get(0,'MonitorPositions');
selected_monpos = monpos(select_monitor,:);
oneplotfig_monpos = selected_monpos;
oneplotfig_monpos(4) = (selected_monpos(4) * 0.8 )/4; % adjust height

for sel_datainfo = all_datainfo
    for resample_fps = all_resample_fps
        for derivative_order = all_derivative_orders

            gran_continue = 1;
            gran_count = 0;
            gran = 0;
            
            for distance_metric_count = 1:length(all_distance_metric) 
                
                while gran_continue
                    
                    distance_metric = all_distance_metric{distance_metric_count};
                    gran_count = gran_count + 1;
                    
                    if isnumeric(all_granularity)
                        if gran_count > length(all_granularity)
                            break
                        else
                            gran = all_granularity(gran_count);
                        end
                    else
                        gran = gran + 1;
                    end
                    
                    for sel_method = all_method
                        if sel_method == 2
                            method_label = 'offline-alternative';
                        elseif sel_method == 1
                            method_label = 'offline-as-online';
                        elseif sel_method == 0
                            method_label = 'online';
                        end
                        di = 0;
                        alldatainfo = struct('skipderiv',[],'gran_continue',gran_continue);
                        
                        % -----------------------------------------------------------------------
                        % INFORMATION ABOUT DATA AND VARIABLES:
                        % novelty_width_s   = free variable novelty kernel (Gaussian-tapered checkerboard) width (in seconds)
                        % filter_width_s    = free variable filter kernel (Gaussian) width (in seconds)
                        % peaks_threshold   = free variable novelty peaks threshold (theta)
                        % segnum_fontsize   = font size for segment numbers
                        % az_el             = azimuth and elevation for pictures (in degrees)
                        % skipderiv         = skip parameters for derivative order (logical)
                        % . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
                        di = di + 1;
                        alldatainfo(di).file_name             = 'man_heaven_six.tsv'; % recording 1
                        alldatainfo(di).xtick_interval        = 2;
                        alldatainfo(di).limits                = [-1.2583e+03, 1.0999e+03, -195.2414, 2.1630e+03];
                        if gran == 1 % coarse granularity: weaving, rest, false starts
                            if derivative_order == 0 % position
                                alldatainfo(di).skipderiv         = 1;      
                            elseif derivative_order == 1 % velocity
                                alldatainfo(di).novelty_width_s   = 4.8;    % <--- 30fps euclidean
                                alldatainfo(di).filter_width_s    = 1.2;    % <---   "       "  
                                alldatainfo(di).peaks_threshold   = 6;      % <---   "       "     discards fine motion
                            end
                            alldatainfo(di).segnum_fontsize   = 30;         % <--- 
                            alldatainfo(di).az_el             = [0, 0];     % <---
                         elseif gran == 2 % fine granularity: left or right hits
                            if derivative_order == 0 % position
                                alldatainfo(di).skipderiv         = 1;      
                            elseif derivative_order == 1 % velocity
                                alldatainfo(di).novelty_width_s   = 2.47;   % <--- 30fps euclidean
                                alldatainfo(di).filter_width_s    = 0.63;   % <---   "       "  
                                alldatainfo(di).peaks_threshold   = 0;      % <---   "       " 
                            end
                            alldatainfo(di).segnum_fontsize   = 10;         % <---
                            alldatainfo(di).az_el             = [180, -90]; % <---
                        else
                            alldatainfo(di).gran_continue = 0;
                        end
                        
                        % . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
                        di = di + 1;
                        alldatainfo(di).file_name             = '14_20.c3d'; % recording 2
                        alldatainfo(di).sample_frequency      = 120;
                        alldatainfo(di).select_markers        = [1, 2, 3, 4, 10, 12, 15, 17, 19, 22, 24, 25, 29, 34, 36, 41];
                        alldatainfo(di).xtick_interval        = 1;
                        alldatainfo(di).limits                = [-6.0173e+22, 6.3311e+22, -8.5356e+21, 1.1495e+23];
                        if gran == 1 % coarse granularity: jump, twist, knees, reach up-down, reach right-left, arm circles
                            if derivative_order == 0 % position
                                alldatainfo(di).skipderiv         = 1;      
                            elseif derivative_order == 1 % velocity
                                alldatainfo(di).novelty_width_s   = 15;     % <--- 30fps euclidean
                                alldatainfo(di).filter_width_s    = 3.7;    % <---   "       "  
                                alldatainfo(di).peaks_threshold   = 0;      % <---   "       " 
                            end
                            alldatainfo(di).segnum_fontsize   = 30;         % <---
                            alldatainfo(di).az_el             = [90,0];     % <---
                        elseif gran == 2 % fine granularity
                            if derivative_order == 0 % position
                                alldatainfo(di).skipderiv         = 1;      
                            elseif derivative_order == 1 % velocity
                                alldatainfo(di).novelty_width_s   = 1.8;    % <--- 30fps euclidean  
                                alldatainfo(di).filter_width_s    = 0.43;   % <---   "       "  
                                alldatainfo(di).peaks_threshold   = 0;      % <---   "       " 
                            end
                            alldatainfo(di).segnum_fontsize   = 10;         % <---
                            alldatainfo(di).az_el             = [90,0];     % <---
                        else
                            alldatainfo(di).gran_continue = 0;
                        end
                        % . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
                        di = di + 1;
                        alldatainfo(di).file_name             = '18_15.c3d'; % recording 3
                        alldatainfo(di).select_markers        = [3:11,13,18,20,22,26,30,32,33,37,41,44,46,47,50:57,60,63,64,67,69,70,72,74,75,80,82]; % subject A: 'Justin'
                        alldatainfo(di).xtick_interval        = 1;
                        alldatainfo(di).limits                = [-1.8599e+03 176.1985 -160.5019 1.8756e+03];
                        if gran == 1 % coarse granularity: beaks, wings, tail, claps
                            if derivative_order == 0 % position     
                                alldatainfo(di).skipderiv         = 1;      
                            elseif derivative_order == 1 % velocity
                                alldatainfo(di).novelty_width_s   = 3.6;    % <--- 30fps euclidean 
                                alldatainfo(di).filter_width_s    = 0.9;    % <---   "       "  
%                                 alldatainfo(di).peaks_threshold   = 3;      % <---   "       "     online or offline
                                alldatainfo(di).peaks_threshold   = 0;      % <---   "       "     only online
                            end
                            alldatainfo(di).segnum_fontsize   = 30;         % <---
                            alldatainfo(di).az_el             = [90, 0];    % <---
                        elseif gran == 2 % fine granularity
                            if derivative_order == 0 % position
                                alldatainfo(di).skipderiv         = 1;      
                            elseif derivative_order == 1 % velocity
                                alldatainfo(di).novelty_width_s   = 0.5;    % <--- 30fps euclidean 
                                alldatainfo(di).filter_width_s    = 0.1;    % <---   "       "  
                                alldatainfo(di).peaks_threshold   = 5;      % <---   "       "     this discards beaks, for lack of markers in fingers
                            end
                            alldatainfo(di).segnum_fontsize   = 10;         % <---
                            alldatainfo(di).az_el             = [90,0];     % <---
                        else
                            alldatainfo(di).gran_continue = 0;
                        end
                        % . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
                        di = di + 1;
                        alldatainfo(di).file_name             = 'man_woman_lindy_hop.tsv'; % recording 4
                        alldatainfo(di).trim                  = [427,15557]; % music stimulus
                        alldatainfo(di).select_markers        = 1:17; % woman
                        alldatainfo(di).xtick_interval        = 10;
                        alldatainfo(di).limits                = [-1.2139e+03 2.2520e+03 -854.3455 2.6116e+03];
                        if gran == 1 % coarse granularity
                            if derivative_order == 0 % position
                                alldatainfo(di).skipderiv         = 1;      
                            elseif derivative_order == 1 % velocity
                                alldatainfo(di).novelty_width_s   = 18;     % <--- 30fps euclidean 
                                alldatainfo(di).filter_width_s    = 4.5;    % <---   "       "  
                                alldatainfo(di).peaks_threshold   = 0;      % <---   "       "    
                            end
                            alldatainfo(di).segnum_fontsize   = 20;         % <---
                            alldatainfo(di).az_el             = [90,0];     % <---
                        elseif gran == 2 % coarse-medium granularity
                            if derivative_order == 0 % position
                                alldatainfo(di).skipderiv         = 1;      
                            elseif derivative_order == 1 % velocity
                                alldatainfo(di).novelty_width_s   = 11;     % <--- 30fps euclidean 
                                alldatainfo(di).filter_width_s    = 2.9;    % <---   "       "  
                                alldatainfo(di).peaks_threshold   = 1.3;    % <---   "       "      
                            end
                            alldatainfo(di).segnum_fontsize   = 20;         % <---
                            alldatainfo(di).az_el             = [90,0];     % <---
                        elseif gran == 3 % medium-fine granularity
                            if derivative_order == 0 % position
                                alldatainfo(di).skipderiv         = 1;      
                            elseif derivative_order == 1 % velocity 
                                alldatainfo(di).novelty_width_s   = 7;      % <--- 30fps euclidean 
                                alldatainfo(di).filter_width_s    = 1.76;   % <---   "       "  
                                alldatainfo(di).peaks_threshold   = 1.9;    % <---   "       "     filters individual kicks
                            end
                            alldatainfo(di).segnum_fontsize   = 10;         % <---
                            alldatainfo(di).az_el             = [90,0];     % <---
                        else
                            alldatainfo(di).gran_continue = 0;
                        end
                        % . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
                        di = di + 1;
                        alldatainfo(di).file_name             = 'woman_disco.tsv'; % recording 5
                        alldatainfo(di).trim                  = [621,11470]; % music stimulus
                        alldatainfo(di).xtick_interval        = 10;
                        alldatainfo(di).limits                = [-1.4278e+03 1.2004e+03 -452.7991 2.1754e+03];
                        if gran == 1 % coarse granularity
                            if derivative_order == 0 % position
                                alldatainfo(di).skipderiv         = 1;      
                            elseif derivative_order == 1 % velocity
                                alldatainfo(di).novelty_width_s   = 8;      % <--- 30fps euclidean
                                alldatainfo(di).filter_width_s    = 2.7;    % <---   "       "       
                                alldatainfo(di).peaks_threshold   = 0;      % <---   "       "     
                            end
                            alldatainfo(di).segnum_fontsize   = 14;         % <---
                            alldatainfo(di).az_el             = [90,0];     % <---
                        elseif gran == 2 % coarse-medium granularity
                            if derivative_order == 0 % position
                                alldatainfo(di).skipderiv         = 1;      
                            elseif derivative_order == 1 % velocity
                                alldatainfo(di).novelty_width_s   = 8;      % <--- 30fps euclidean
                                alldatainfo(di).filter_width_s    = alldatainfo(di).novelty_width_s / 4;    % <---   "       "    
                                alldatainfo(di).peaks_threshold   = 0;      % <---   "       "     
                            end
                            alldatainfo(di).segnum_fontsize   = 10;         % <---
                            alldatainfo(di).az_el             = [90,0];     % <---
                        elseif gran == 3 % medium-fine granularity
                            if derivative_order == 0 % position
                                alldatainfo(di).skipderiv         = 1;      
                            elseif derivative_order == 1 % velocity
                                alldatainfo(di).novelty_width_s   = 6.6;    % <--- 30fps euclidean
                                alldatainfo(di).filter_width_s    = alldatainfo(di).novelty_width_s / 4;    % <---   "       "    
                                alldatainfo(di).peaks_threshold   = 0;      % <---   "       "     
                            end
                            alldatainfo(di).segnum_fontsize   = 10;         % <---
                            alldatainfo(di).az_el             = [90,0];     % <---
                        elseif gran == 4 % fine granularity
                            if derivative_order == 0 % position
                                alldatainfo(di).skipderiv         = 1;      
                            elseif derivative_order == 1 % velocity
                                alldatainfo(di).novelty_width_s   = 6.6;    % <--- 30fps euclidean
                                alldatainfo(di).filter_width_s    = 1.1;    % <---   "       "   
                                alldatainfo(di).peaks_threshold   = 0;      % <---   "       "     
                            end
                            alldatainfo(di).segnum_fontsize   = 10;         % <---
                            alldatainfo(di).az_el             = [90,0];     % <---
                        else
                            alldatainfo(di).gran_continue = 0;
                        end
                        % . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
                        % add more recordings using the same logic
                        % -----------------------------------------------------------------------
                        
                        clear datainfo
                        datainfo = alldatainfo(sel_datainfo);
                        
                        % Make derivative labels:
                        if derivative_order == 0
                            derivative_label = 'position';
                        elseif derivative_order == 1
                            derivative_label = 'velocity';
                        elseif derivative_order == 2
                            derivative_label = 'acceleration';
                        else
                            error('Only derivative orders 0, 1 and 2 are allowed')
                        end

                        if isempty(datainfo.skipderiv) == 0 
                            disp(sprintf('skipping segmentation for recording %g at granularity level = %d and derivative order = %d',sel_datainfo,gran,derivative_order))
                        else
                            if datainfo.gran_continue == 0
                                gran_continue = 0;
                                break
                            end
                            if sel_datainfo > di
                                error('At least one selected recording is out of range.')
                            end
                            
                            if isempty(economy_sw) == 0
                                
                                % To make computation faster, the structures with "_change" in their names
                                % track changes in variables for critical steps, so those steps are not
                                % done if their results have already been computed.
                                if exist('import_change','var') == 0 ...
                                        || strcmp(import_change.file_name,datainfo.file_name) == 0
                                    
                                    % Import mocap data:
                                    cd(data_path)
                                    mocap_struct_raw = mcread(datainfo.file_name );
                                    if isfield(datainfo,'sample_frequency') && isempty(datainfo.sample_frequency) == 0
                                        mocap_struct_raw.freq = datainfo.sample_frequency;
                                    end
                                    [~,file_name,~] = fileparts(datainfo.file_name);
                                    animpar = mcinitanimpar;
                                    
                                    import_change.file_name = datainfo.file_name;
                                end
                                
                                if isfield(datainfo,'trim') == 0
                                    datainfo.trim = [];
                                end
                                
                                if isempty(datainfo.trim)
                                    mocap_struct_trimmed = mocap_struct_raw;
                                elseif (exist('trim_change','var') == 0 ...
                                        || strcmp(trim_change.file_name,datainfo.file_name) == 0 ...
                                        || trim_change.resample_fps ~= resample_fps ...
                                        || isequal(trim_change.trim, datainfo.trim) == 0)
                                    
                                    % Trim  mocap data:
                                    cd(data_path)
                                    mocap_struct_trimmed = mctrim(mocap_struct_raw, datainfo.trim(1), datainfo.trim(2), 'frame' );
                                    
                                    trim_change.resample_fps = resample_fps;
                                    trim_change.file_name = datainfo.file_name;
                                    trim_change.trim = datainfo.trim;
                                end

                                if exist('preprocess_change','var') == 0 ...
                                        || strcmp(preprocess_change.file_name,datainfo.file_name) == 0 ...
                                        || preprocess_change.resample_fps ~= resample_fps ...
                                        || isequal(preprocess_change.trim, datainfo.trim) == 0 ...
                                        || isequal(preprocess_change.select_markers,datainfo.select_markers) == 0 ...
                                        || isequal(preprocess_change.derivative_order,derivative_order) == 0
                                    
                                    % Extract selected markers:
                                    if isfield(datainfo,'select_markers') && isempty(datainfo.select_markers) == 0
                                        mocap_struct_selected = mcgetmarker(mocap_struct_trimmed, datainfo.select_markers);
                                    else
                                        mocap_struct_selected = mocap_struct_trimmed;
                                    end
                                    
                                    % Fill gaps:
                                    mocap_struct_selected = mcfillgaps(mocap_struct_selected);
                                    
                                    if isempty(all_selected_segments) == 0

                                        % Resample:
                                        if mocap_struct_selected.freq ~= resample_fps
                                            mocap_struct_resampled = mcresample(mocap_struct_selected,resample_fps);
                                        else
                                            mocap_struct_resampled = mocap_struct_selected;
                                        end
                                    end
                                    
                                    preprocess_change.file_name = datainfo.file_name;
                                    preprocess_change.trim = datainfo.trim;
                                    preprocess_change.select_markers = datainfo.select_markers;
                                    preprocess_change.derivative_order = derivative_order;
                                    preprocess_change.resample_fps = resample_fps;
                                end
                            end
                            
                            if isempty(all_selected_segments) == 0
                                
                                datainfo.novelty_width_f = round(datainfo.novelty_width_s * resample_fps); % novelty kernel width in frames = n(nov)
                                if rem(datainfo.novelty_width_f,2) % ensure that the width of the novelty kernel is even so that the checkerboard can be produced
                                    datainfo.novelty_width_f = datainfo.novelty_width_f + 1;
                                end
                                
                                datainfo.filter_width_f = round(datainfo.filter_width_s * resample_fps); % filter kernel width in frames = n(filt)
                                if rem(datainfo.filter_width_f,2) == 0 % ensure that the width of the filter kernel is odd so that gaussian is centered
                                    datainfo.filter_width_f = datainfo.filter_width_f - 1;
                                end
                                if datainfo.filter_width_f <= 1
                                    datainfo.filter_width_f = 0;
                                end
 
                                distance_label = distance_metric(1:3);

                                if isempty(economy_sw) == 0

                                    % ------------------------------------------------------------------
                                    % KERNELS
                                    
                                    % 2D Gaussian-tapered checkerboard kernel for novelty:
                                    if exist('novelty_params_change','var') == 0 ...
                                            || novelty_params_change.novelty_width_f ~= datainfo.novelty_width_f ...
                                            || economy_sw == 0
                                        
                                        novelty_alpha = 2.5; % <--- times the sdev fits into half the novelty kernel's width
                                        
                                        [novelty_kernel_x, novelty_kernel_y] = meshgrid((-(datainfo.novelty_width_f-1)/2):((datainfo.novelty_width_f-1)/2), (-(datainfo.novelty_width_f-1)/2):((datainfo.novelty_width_f-1)/2));
                                        novelty_kernel.variance = ((datainfo.novelty_width_f - 1 )/(2*novelty_alpha)) ^2;
                                        novelty_kernel.gauss = exp( - ( (novelty_kernel_x.^2) / ( 2 * novelty_kernel.variance) ) - ( (novelty_kernel_y.^2) / ( 2 * novelty_kernel.variance) ) );
                                        novelty_kernel.gauss = novelty_kernel.gauss - min(novelty_kernel.gauss(:)); % minima zero
                                        novelty_kernel.gauss = novelty_kernel.gauss / sum(novelty_kernel.gauss(:)); % unit volume
                                        novelty_kernel.cb = (kron([-1, 1; 1,  -1],ones(datainfo.novelty_width_f/2))); % checkerboard
                                        novelty_kernel.gausscb = novelty_kernel.cb .* novelty_kernel.gauss; % gaussian-tapered checkerboard
                                        
                                        novelty_params_change.novelty_width_f = datainfo.novelty_width_f;
                                    end
                                    
                                    % 1D Gaussian kernel for filtering:
                                    if datainfo.filter_width_f > 1 ...
                                            && exist('filter_params_change','var') == 0 ...
                                            || filter_params_change.filter_width_f ~= datainfo.filter_width_f ...
                                            || economy_sw == 0
                                        
                                        filter_alpha = 2.5; % <--- times the sdev fits into half the filter kernel's width
                                        
                                        filter_kernel_x = linspace( -(datainfo.filter_width_f-1)/2, (datainfo.filter_width_f-1)/2, datainfo.filter_width_f )';
                                        filter_kernel.variance = ((datainfo.filter_width_f - 1 )/(2*filter_alpha))^2;
                                        filter_kernel.gauss = exp( -(filter_kernel_x.^2) / (2*filter_kernel.variance) );
                                        filter_kernel.gauss = filter_kernel.gauss - min(filter_kernel.gauss(:)); % confine
                                        filter_kernel.gauss = filter_kernel.gauss / sum(filter_kernel.gauss(:)); % unit area
                                        
                                        filter_params_change.filter_width_f = datainfo.filter_width_f;
                                    end
                                    
                                    if sel_method == 0 || sel_method == 1
                                        lag_offset = datainfo.novelty_width_f/2 + ceil(datainfo.filter_width_f/2) + 1; % <--- offset to compensate lags of novelty, filter and differential (in frames)
                                    elseif sel_method == 2
                                        if datainfo.filter_width_f > 1
                                            lag_offset = ceil(datainfo.filter_width_f/2) + 1; % <--- filter applied: offset to compensate lags of filter and differential (in frames)
                                        else
                                            lag_offset =  1; % <--- filter not applied: offset to compensate differential (in frames)
                                        end
                                    end
                                    
                                    process_track = 0;
                                    if sel_method == 0
                                        
                                        % ------------------------------------------------------------------
                                        % ON-LINE SEGMENTATION
                                        
                                        if exist('online_seg_change','var') == 0 ...
                                                || strcmp(online_seg_change.file_name,datainfo.file_name) == 0 ...
                                                || isequal(online_seg_change.trim, datainfo.trim) == 0 ...
                                                || isequal(online_seg_change.select_markers,datainfo.select_markers) == 0 ...
                                                || online_seg_change.resample_fps ~= resample_fps ...
                                                || isequal(online_seg_change.derivative_order,derivative_order) == 0 ...
                                                || strcmp(online_seg_change.distance_metric, distance_metric) == 0 ...
                                                || online_seg_change.novelty_width_f ~= datainfo.novelty_width_f ...
                                                || online_seg_change.filter_width_f ~= datainfo.filter_width_f ...
                                                || online_seg_change.peaks_threshold ~= datainfo.peaks_threshold ...
                                                || exist('sel_method_change','var') == 0 || sel_method_change ~= sel_method ...
                                                || economy_sw == 0
                                            
                                            disp([method_label,' computation'])
                                            process_track = 1;
                                            
                                            % make upper triangle of novelty kernel into a vector:
                                            novelty_kernel.gausscb_triu = triu(novelty_kernel.gausscb,1);
                                            novelty_kernel.gausscb_triu_i = triu(true(datainfo.novelty_width_f,datainfo.novelty_width_f),1);
                                            novelty_kernel.gausscb_triu_v = flipud(novelty_kernel.gausscb_triu(novelty_kernel.gausscb_triu_i(:)));
                                            
                                            % Initialise buffers:
                                            if derivative_order > 0
                                                differentiation_buffer = zeros(derivative_order+1,size(mocap_struct_resampled.data,2));
                                            end
                                            frames_buffer_length = datainfo.novelty_width_f * 8; % <--- length for frames buffer (novelty is computed with values in the novelty window within this buffer)
                                            novelty_buffer_length = datainfo.filter_width_f * 8; % <--- length for novelty buffer (filtered novelty is computed with values in the filter window within this buffer)
                                            
                                            frames_buffer = zeros(frames_buffer_length,mocap_struct_resampled.nMarkers*3);
                                            novelty_buffer = zeros(1,novelty_buffer_length);
                                            if strcmp(distance_metric,'euclidean')
                                                online_distmat_triu = zeros(datainfo.novelty_width_f); % distance matrix buffer
                                                online_distmat_triu_i = novelty_kernel.gausscb_triu_i; % index for distance matrix's upper triangle
                                                online_distmat_triu_v = zeros(1,(datainfo.novelty_width_f^2 - datainfo.novelty_width_f)/2); % initialise vector of distance matrix's upper triangle
                                            end
                                            peaktest_buffer = [0,0,0];
                                            
                                            % Initialise other variables:
                                            nov_filt = NaN(1, mocap_struct_resampled.nFrames + lag_offset); % filtered novelty score
                                            seg_bounds_count = 0; % counter for novelty peaks (segmentation boundaries)
                                            seg_bounds = nov_filt; % vector to store the maximum possible amount of bondaries (data length)
                                            
                                            % Initialise counters:
                                            novelty_window_begin = 0;
                                            novelty_window_end = 0;
                                            frames_buffer_margin_index = 0;
                                            filter_window_begin = 0;
                                            filter_window_end = 0;
                                            novelty_buffer_margin_index = 0;
                                            
                                            % initialise figure for real-time display of online process:
                                            if plotprocess_sw
                                                
                                                rtfig_monpos = selected_monpos;
                                                rtfig_monpos(4) = rtfig_monpos(4)/2;
                                                figh_olprocplot = figure('Position',rtfig_monpos);
                                                rt_title_fontsize      = 22;
                                                rt_ticklabels_fontsize = 14;
                                                rt_axeslabels_fontsize = 16;
                                                rt_boundaries_colour = [0.2 0.7 0.2];
                                                rt_threshold_colour = ([1 1 1] * 0.6) - 0.5;
                                                rt_displayed_seg_bounds_count = 0;
                                                rt_new_YTick = [0,0,0];
                                                rt_max_nov_filt = 0;
                                                
                                                % mocap data buffer plot:
                                                subplot(2,3,1)
                                                ph_mcdata_buffer = plot(sum(frames_buffer,2));
                                                set(gca,'fontsize',rt_ticklabels_fontsize)
                                                xlabel('frames','fontsize',rt_axeslabels_fontsize)
                                                ylabel('sum','fontsize',rt_axeslabels_fontsize)
                                                title('mocap data buffer','fontsize',rt_title_fontsize,'fontweight','normal')
                                                
                                                % self-similarity matrix image:
                                                subplot(2,3,2)
                                                ph_distmat = imagesc(squareform(online_distmat_triu_v));
                                                set(gca,'fontsize',rt_ticklabels_fontsize)
                                                xlabel('frames','fontsize',rt_axeslabels_fontsize)
                                                ylabel('frames','fontsize',rt_axeslabels_fontsize)
                                                axis square
                                                title('online self-similarity matrix','fontsize',rt_title_fontsize,'fontweight','normal')
                                                
                                                % novelty buffer plot:
                                                subplot(2,3,3)
                                                ph_novelty_buffer = plot(novelty_buffer);
                                                set(gca,'fontsize',rt_ticklabels_fontsize)
                                                xlabel('frames','fontsize',rt_axeslabels_fontsize)
                                                ylabel('novelty','fontsize',rt_axeslabels_fontsize)
                                                title('novelty buffer','fontsize',rt_title_fontsize,'fontweight','normal')
                                                
                                                % filtered novelty  plot:
                                                subplot(2,1,2)
                                                ph_nov_filt = plot(nov_filt);
                                                ph_nov_filt.Parent.Position(4) = ph_nov_filt.Parent.Position(4) * 0.8;
                                                rt_xticks = [0:datainfo.xtick_interval*resample_fps:length(nov_filt)+lag_offset];
                                                rt_xticklabels = [0 : datainfo.xtick_interval : (length(nov_filt)+lag_offset)/resample_fps ];
                                                set(gca,'xlim',[1,length(nov_filt)],'xtick',rt_xticks,'xticklabel',rt_xticklabels,...
                                                    'Ytick',unique(sort([0,datainfo.peaks_threshold])),'fontsize',rt_ticklabels_fontsize)
                                                xlabel('time (seconds)','fontsize',rt_axeslabels_fontsize)
                                                ylabel('novelty','fontsize',rt_axeslabels_fontsize)
                                                title('filtered novelty score','fontsize',rt_title_fontsize,'fontweight','normal')
                                                
                                                % plot peak threshold line:
                                                hold on
                                                pkthr_line = ones(1,length(nov_filt)+lag_offset) * datainfo.peaks_threshold;
                                                plot(pkthr_line,'--','color',rt_threshold_colour,'LineWidth',1)
                                                hold off
                                            end
                                            
                                            % start online process %
                                            % note that the system of buffers is intended for processing one mocap frame at a time, which suits real-time application
                                            process_start_time = tic;
                                            for online_seg_count = 1:mocap_struct_resampled.nFrames + lag_offset
                                                
                                                % Compute beginning and ending indexes of current window within the mocap data buffer:
                                                if novelty_window_end == frames_buffer_length
                                                    novelty_window_begin = 0;
                                                    frames_buffer_margin_index = 0;
                                                end
                                                novelty_window_begin = novelty_window_begin + 1;
                                                novelty_window_end = novelty_window_begin + datainfo.novelty_width_f - 1;

                                                % Populate data buffer:
                                                if online_seg_count <= mocap_struct_resampled.nFrames % stop at end of mocap data
                                                    incoming_frame = mocap_struct_resampled.data( online_seg_count , : );
                                                    if derivative_order > 0 % Derivate by difference:
                                                        differentiation_buffer(1:derivative_order,:) = differentiation_buffer(2:end,:); % shift differentiation buffer
                                                        differentiation_buffer(end,:) = incoming_frame; % add mocap data frame at the end of the buffer
                                                        if online_seg_count <= derivative_order
                                                            incoming_frame = zeros(1,mocap_struct_resampled.nMarkers * 3); % omit differences with zeroes in buffer 
                                                        else
                                                            incoming_frame = diff(differentiation_buffer,derivative_order);
                                                        end
                                                    end
                                                    frames_buffer(novelty_window_end,:) = incoming_frame; % add new data vector at the end of the buffer
                                                    if (frames_buffer_length - datainfo.novelty_width_f +1) < novelty_window_end % populate the beginning after first reset
                                                        frames_buffer_margin_index = frames_buffer_margin_index + 1;
                                                        frames_buffer(frames_buffer_margin_index,:) = frames_buffer(novelty_window_end,:);
                                                    end
                                                else % after the end of mocap data populate with zeros
                                                    frames_buffer(novelty_window_end,:) = zeros(1,mocap_struct_resampled.nMarkers * 3);
                                                end
                                                
                                                % Distance matrix of frames buffer:
                                                if strcmp(distance_metric,'euclidean')
                                                    online_distmat_triu(2:end,2:end) = online_distmat_triu(1:end-1,1:end-1); % shift
                                                    for euclidean_counter = 1:datainfo.novelty_width_f
                                                        online_distmat_triu(1,euclidean_counter) = ...
                                                            sqrt(sum((frames_buffer(novelty_window_end,:) - frames_buffer(novelty_window_end-euclidean_counter+1,:)) .^2));
                                                    end
                                                    online_distmat_triu_v = flipud(online_distmat_triu(online_distmat_triu_i(:)))';
                                                else
                                                    online_distmat_triu_v = pdist(frames_buffer(novelty_window_begin:novelty_window_end,:),distance_metric);
                                                end
                                                
                                                % Compute beginning and ending indexes of current window within the novelty buffer:
                                                if filter_window_end == novelty_buffer_length
                                                    filter_window_begin = 0;
                                                    novelty_buffer_margin_index = 0;
                                                end
                                                filter_window_begin = filter_window_begin + 1;
                                                filter_window_end = filter_window_begin + datainfo.filter_width_f - 1;
                                                
                                                % Compute novelty and put result at the end of the window in the novelty buffer:
                                                novelty_buffer(filter_window_end) = online_distmat_triu_v * novelty_kernel.gausscb_triu_v; % scalar product of distance matrix with novelty kernel
                                                if (novelty_buffer_length - datainfo.filter_width_f + 1) < filter_window_end % populate the beginning after first reset
                                                    novelty_buffer_margin_index = novelty_buffer_margin_index + 1;
                                                    novelty_buffer(novelty_buffer_margin_index) = novelty_buffer(filter_window_end);
                                                end
                                                
                                                % Shift peak test buffer:
                                                peaktest_buffer(1:2) = peaktest_buffer(2:3);
                                                
                                                % Apply Gaussian filter to the filter buffer and put result at the end of the difference buffer:
                                                if datainfo.filter_width_f > 1
                                                    peaktest_buffer(end) = novelty_buffer(filter_window_begin:filter_window_end) * filter_kernel.gauss; % scalar product of filter buffer with half-gauss filter kernel
                                                else
                                                    peaktest_buffer(end) = novelty_buffer(filter_window_end); % don't filter if width of filter kernel is one frame
                                                end
                                                nov_filt(online_seg_count) = peaktest_buffer(end);
                                                
                                                % Extract peaks over threshold:
                                                [~,local_maxima_test] = max(peaktest_buffer); 
                                                if local_maxima_test == 2 && peaktest_buffer(2) >= datainfo.peaks_threshold 
                                                        seg_bounds_count = seg_bounds_count + 1;
                                                        seg_bounds(seg_bounds_count) = online_seg_count; % main output
                                                end                                                
                                                
                                                % Plot process in real time:
                                                if plotprocess_sw
                                                    set(ph_mcdata_buffer,'YData',(sum(frames_buffer,2)))
                                                    set(ph_distmat,'CData',squareform(online_distmat_triu_v))
                                                    set(ph_novelty_buffer,'YData',novelty_buffer)
                                                    set(ph_nov_filt,'YData',nov_filt)
                                                    if peaktest_buffer(end) > rt_max_nov_filt
                                                        rt_max_nov_filt = peaktest_buffer(end);
                                                        rt_new_YTick = unique(sort([min(nov_filt),datainfo.peaks_threshold,rt_max_nov_filt]));
                                                        set(ph_nov_filt.Parent,'YTick',rt_new_YTick)
                                                        if rt_new_YTick(end) > datainfo.peaks_threshold
                                                            set(ph_nov_filt.Parent,'YLim',[0,rt_new_YTick(end)])
                                                        end
                                                    end
                                                    if seg_bounds_count > rt_displayed_seg_bounds_count
                                                        rt_displayed_seg_bounds_count = seg_bounds_count;
                                                        rep_peaks = repmat(seg_bounds(seg_bounds>0)-1,2,1);
                                                        lines = repmat(get(gca,'ylim'),size(rep_peaks,2),1)';
                                                        line(rep_peaks,lines,'Color',rt_boundaries_colour,'linewidth',3,'linestyle','-');
                                                    end
                                                    drawnow
                                                end                                                
                                            end
                                            process_stop_time = toc(process_start_time);
                                            % end online process %
                                            
                                            seg_bounds = seg_bounds(1:seg_bounds_count); % trim to keep only found boundaries
                                            if lag_offset % compensate lag for assessment
                                                nov_filt = nov_filt(lag_offset:lag_offset+mocap_struct_resampled.nFrames-1);
                                                seg_bounds = seg_bounds - lag_offset;
                                            end
                                            seg_bounds = seg_bounds(seg_bounds > 0); % eliminate any boundaries that may be slightly before zero time
                                            
                                            online_seg_change.file_name = datainfo.file_name;
                                            online_seg_change.trim = datainfo.trim;
                                            online_seg_change.select_markers = datainfo.select_markers;
                                            online_seg_change.resample_fps = resample_fps;
                                            online_seg_change.derivative_order = derivative_order;
                                            online_seg_change.distance_metric = distance_metric;
                                            online_seg_change.novelty_width_f = datainfo.novelty_width_f;
                                            online_seg_change.filter_width_f = datainfo.filter_width_f;
                                            online_seg_change.peaks_threshold = datainfo.peaks_threshold;
                                        end
                                    elseif sel_method == 1 || sel_method == 2
                                        
                                        % ------------------------------------------------------------------
                                        % OFF-LINE SEGMENTATION
                                        
                                        if exist('offline_nov_change','var') == 0 ...
                                                || strcmp(offline_nov_change.file_name,datainfo.file_name) == 0 ...
                                                || isequal(offline_nov_change.trim, datainfo.trim) == 0 ...
                                                || isequal(offline_nov_change.select_markers,datainfo.select_markers) == 0 ...
                                                || offline_nov_change.resample_fps ~= resample_fps ...
                                                || isequal(offline_nov_change.derivative_order, derivative_order) == 0 ...
                                                || strcmp(offline_nov_change.distance_metric, distance_metric) == 0 ...
                                                || offline_nov_change.novelty_width_f ~= datainfo.novelty_width_f...
                                                || offline_nov_change.filter_width_f ~= datainfo.filter_width_f...
                                                || exist('sel_method_change','var') == 0 || sel_method_change ~= sel_method...
                                                || economy_sw == 0
                                            
                                            disp([method_label,' computation of novelty'])
                                            process_track = 2;
                                            
                                            % Derivate by difference:
                                            if derivative_order > 0
                                                mocap_struct_differentiated = mocap_struct_resampled;
                                                mocap_struct_differentiated.data(derivative_order+1:end-derivative_order+1,:) = diff(mocap_struct_differentiated.data,derivative_order);
                                                mocap_struct_differentiated.data(1:derivative_order,:) = zeros(1,mocap_struct_resampled.nMarkers * 3); % pad with zeroes (emulate "omit differences with zeroes in buffer")
                                                mocap_struct_differentiated.timederOrder = derivative_order;
                                            end
                                            
                                            % initialise vectors:
                                            nov = zeros(1, mocap_struct_differentiated.nFrames + 2); % add 2 samples at ending for testing local maxima
                                            
                                            % start offline process %
                                            process_start_time = tic;
                                            
                                            if sel_method == 1 % offline as online
                                                
                                                % Pad beginning and ending of data:
                                                extended_mcdata = zeros(datainfo.novelty_width_f*2 + mocap_struct_differentiated.nFrames,...
                                                    mocap_struct_differentiated.nMarkers * 3);
                                                extended_mcdata(datainfo.novelty_width_f + 1:end-datainfo.novelty_width_f,:) = mocap_struct_differentiated.data;
                                                
                                                % Distance matrix:
                                                offline_distmat = squareform(pdist(extended_mcdata,distance_metric));
                                                
                                                % Correlation length:
                                                corr_length = mocap_struct_differentiated.nFrames + lag_offset;
                                                
                                            elseif sel_method == 2 % offline alternative
                                                
                                                % Distance matrix:
                                                offline_distmat_raw = mcsimmat(mocap_struct_differentiated,distance_metric);
                                                
                                                offline_beginning = datainfo.novelty_width_f/2 + 1;
                                                offline_ending = datainfo.novelty_width_f/2 + mocap_struct_differentiated.nFrames;
                                                offline_distmat = zeros(mocap_struct_differentiated.nFrames + datainfo.novelty_width_f);
                                                
                                                % Pad beginning and ending corners of distance matrix with averages:
                                                offline_distmat(1:datainfo.novelty_width_f,1:datainfo.novelty_width_f) = ...
                                                    mean(mean(offline_distmat_raw(1:datainfo.novelty_width_f/2,1:datainfo.novelty_width_f/2)));
                                                offline_distmat(end-datainfo.novelty_width_f:end,end-datainfo.novelty_width_f:end) =...
                                                    mean(mean(offline_distmat_raw(end-datainfo.novelty_width_f/2:end,end-datainfo.novelty_width_f/2:end)));
                                                offline_distmat( offline_beginning : offline_ending , offline_beginning : offline_ending ) = ...
                                                    offline_distmat_raw;
                                                
                                                % Correlation length:
                                                corr_length = mocap_struct_differentiated.nFrames + 1;
                                            end
                                            
                                            % Correlate checkerboard kernel along diagonal of distance matrix:
                                            for offline_window_begin = 1:corr_length
                                                offline_window_end = offline_window_begin + datainfo.novelty_width_f - 1 ;
                                                this_offline_distmat = offline_distmat( offline_window_begin : offline_window_end , offline_window_begin : offline_window_end );
                                                nov(offline_window_begin) = this_offline_distmat(:)' * novelty_kernel.gausscb(:) ;
                                            end
                                            
                                            nov = nov / 2; % scale to the online computation
                                            nov(nov < 0) = 0; % get rid of negative novelty values
                                            
                                            % Apply Gaussian filter:
                                            if datainfo.filter_width_f > 1
                                                nov_filt = conv(nov,filter_kernel.gauss,'full');
                                            else
                                                nov_filt = nov; % don't filter if width of filter kernel is one frame
                                            end
                                            
                                            % trim and compensate delay:
                                            nov_filt_plus = nov_filt(lag_offset+1:end);
                                            nov_filt = nov_filt_plus(1:mocap_struct_differentiated.nFrames); % length of resampled data
                                            
                                            offline_nov_change.file_name = datainfo.file_name;
                                            offline_nov_change.trim = datainfo.trim;
                                            offline_nov_change.select_markers = datainfo.select_markers;
                                            offline_nov_change.resample_fps = resample_fps;
                                            offline_nov_change.derivative_order = derivative_order;
                                            offline_nov_change.distance_metric = distance_metric;
                                            offline_nov_change.novelty_width_f = datainfo.novelty_width_f;
                                            offline_nov_change.filter_width_f = datainfo.filter_width_f;
                                        end
                                        
                                        if exist('offline_pthr_change','var') == 0 ...
                                                || strcmp(offline_pthr_change.file_name,datainfo.file_name) == 0 ...
                                                || isequal(offline_pthr_change.trim, datainfo.trim) == 0 ...
                                                || isequal(offline_pthr_change.select_markers,datainfo.select_markers) == 0 ...
                                                || offline_pthr_change.resample_fps ~= resample_fps ...
                                                || isequal(offline_pthr_change.derivative_order,derivative_order) == 0 ...
                                                || strcmp(offline_pthr_change.distance_metric, distance_metric) == 0 ...
                                                || offline_pthr_change.novelty_width_f ~= datainfo.novelty_width_f ...
                                                || offline_pthr_change.filter_width_f ~= datainfo.filter_width_f ...
                                                || offline_pthr_change.peaks_threshold ~= datainfo.peaks_threshold ...
                                                || exist('sel_method_change','var') == 0 || sel_method_change ~= sel_method ...
                                                || economy_sw == 0
                                            
                                            disp([method_label,' computation of peaks'])
                                            process_track = 3;
                                            
                                            % Extract peaks over threshold:
                                            peaks_index_all = zeros(1,length(nov_filt_plus));
                                            peaks_index_all(2:end-1) = (diff(sign(diff(nov_filt_plus))) == -2);
                                            peaks_values_all = nov_filt_plus .* peaks_index_all;
                                            peaks_values_all(peaks_values_all == 0) = NaN;
                                            peaks_values_selected = (peaks_values_all >= datainfo.peaks_threshold );
                                            seg_bounds = find( peaks_values_selected == 1 );
                                            
                                            offline_pthr_change.file_name = datainfo.file_name;
                                            offline_pthr_change.trim = datainfo.trim;
                                            offline_pthr_change.select_markers = datainfo.select_markers;
                                            offline_pthr_change.resample_fps = resample_fps;
                                            offline_pthr_change.derivative_order = derivative_order;
                                            offline_pthr_change.distance_metric = distance_metric;
                                            offline_pthr_change.novelty_width_f = datainfo.novelty_width_f;
                                            offline_pthr_change.filter_width_f = datainfo.filter_width_f;
                                            offline_pthr_change.peaks_threshold = datainfo.peaks_threshold;
                                        end
                                        process_stop_time = toc(process_start_time);
                                        % end offline process %
                                        
                                        if plotprocess_sw
                                            % plot self-similarity matrix:
                                            imagesc(offline_distmat)
                                            ssm_title_fontsize      = 22;
                                            ssm_ticklabels_fontsize = 16;
                                            ssm_axeslabels_fontsize = 16;
                                            ssm_xticks = [0:datainfo.xtick_interval*resample_fps:mocap_struct_differentiated.nFrames+datainfo.novelty_width_f*2];
                                            ssm_xticklabels = ...
                                                round([-datainfo.novelty_width_f/resample_fps  : datainfo.xtick_interval : (mocap_struct_differentiated.nFrames+datainfo.novelty_width_f*2)/resample_fps ]);
                                            set(gca,...
                                                'yticklabel',[],...
                                                'xtick',ssm_xticks,...
                                                'xticklabel',ssm_xticklabels)
                                            xlabel('time (seconds)','fontsize',ssm_axeslabels_fontsize)
                                            title([method_label,' self-similarity matrix'],'fontsize',ssm_title_fontsize,'fontweight','normal')
                                        end
                                    end
                                    
                                    % Display computation times or usage of values in memory:
                                    if process_track == 0
                                        disp(sprintf('No segmentation computed. Using %s-computed values in memory.',method_label));
                                    else
                                        timer_minutes = floor(process_stop_time / 60);
                                        timer_seconds = process_stop_time - timer_minutes * 60;
                                        disp(sprintf('process time: %2d minutes, %g seconds', timer_minutes, timer_seconds))
                                    end
                                    
                                    results_label_seg = sprintf('%s_%s_nov%g_filt%g_pt%g',...
                                        derivative_label(1:3),distance_label,datainfo.novelty_width_s,datainfo.filter_width_s,datainfo.peaks_threshold);
                                    
                                    sel_method_change = sel_method;
                                    
                                    if segplot_sw || plotsaveallresults_sw
                                        plotinfo.min_nov_filt(results_count+1) = min(nov_filt);
                                        plotinfo.max_nov_filt(results_count+1) = max(nov_filt);
                                        plotinfo.xticks{results_count+1}       = [0:datainfo.xtick_interval*resample_fps:mocap_struct_resampled.nFrames];
                                        plotinfo.xticklabels{results_count+1}  = [0 : datainfo.xtick_interval : mocap_struct_resampled.nFrames/resample_fps ];
                                        plotinfo.yticks{results_count+1} = [ plotinfo.min_nov_filt(results_count+1) , 0, plotinfo.max_nov_filt(results_count+1) ];
                                        plotinfo.yticks{results_count+1} = sort(plotinfo.yticks{results_count+1});
                                        plotinfo.yticks{results_count+1} = unique(plotinfo.yticks{results_count+1});
                                    end
                                    
                                    % ----------------------------------------------------------------------
                                    % SEGMENTATION PLOT FOR ONE SEQUENCE
                                    
                                    if segplot_sw
                                        if exist('preview_sw','var') == 0 || ( exist('preview_sw','var') && isempty(preview_frame_time) )
                                           
                                            seg_bounds_s = seg_bounds/resample_fps;
                                            disp(sprintf('\n%g SEGMENTATION BOUNDARIES (indexes in seconds):',length(seg_bounds_s)))
                                            disp(sprintf('%g  ',round(seg_bounds_s,2)))
                                            disp(' ')
                                            
                                            if dslabels_sw
                                                dslabel = sprintf('recording %g',sel_datainfo);
                                            else
                                                dslabel = file_name;
                                            end
                                            
                                            figh_segplot = figure('Position',oneplotfig_monpos);
                                            
                                            % plot filtered novelty score:
                                            plot(nov_filt,'color',seg_curve_colour)
                                            set(gca,...
                                                'xlim',[1,mocap_struct_resampled.nFrames],...
                                                'xtick',plotinfo.xticks{results_count+1},...
                                                'xticklabel',plotinfo.xticklabels{results_count+1},...
                                                'ylim',[ plotinfo.min_nov_filt(results_count+1) , plotinfo.max_nov_filt(results_count+1) ],...
                                                'ytick',plotinfo.yticks{results_count+1},...
                                                'yticklabel',sprintf('%1.2g\n',plotinfo.yticks{results_count+1}),...
                                                'fontsize',oneplot_ticklabels_fontsize)
                                            xlabel('time (seconds)','fontsize',oneplot_axeslabels_fontsize)
                                            ylabel('novelty','fontsize',oneplot_axeslabels_fontsize)
                                            hold on
                                            
                                            % plot computed boundaries:
                                            rep_peaks = repmat(seg_bounds,2,1);
                                            lines = repmat(get(gca,'ylim'),size(rep_peaks,2),1)';
                                            line(rep_peaks,lines,'Color',seg_boundaries_colour,'linewidth',3,'linestyle','-');
                                            
                                            % make title:
                                            if rem( datainfo.novelty_width_f, datainfo.novelty_width_s)
                                                approx_nov_ch = sprintf('\\approx');
                                            else
                                                approx_nov_ch = '';
                                            end
                                            if datainfo.novelty_width_f > 1
                                                nov_frames_plural = 's';
                                            else
                                                nov_frames_plural = '';
                                            end
                                            if rem( datainfo.filter_width_f, datainfo.filter_width_s) == 0
                                                approx_filt_ch = '';
                                            else
                                                approx_filt_ch = sprintf('\\approx');
                                            end
                                            if datainfo.filter_width_f > 1
                                                filt_frames_plural = 's';
                                            else
                                                filt_frames_plural = '';
                                            end
                                            title(sprintf('%s , %g fps , %s , %s distance , \\it n\\rm_{nov} = %g frame%s (%s%g s) , \\it n\\rm_{filt} = %g frame%s (%s%g s) , \\theta = %g , %s',...
                                                dslabel, mocap_struct_resampled.freq, derivative_label, distance_metric, ...
                                                datainfo.novelty_width_f, nov_frames_plural, approx_nov_ch, datainfo.novelty_width_s, ...
                                                datainfo.filter_width_f, filt_frames_plural, approx_filt_ch, datainfo.filter_width_s, ...
                                                datainfo.peaks_threshold,method_label),...
                                                'fontsize',oneplot_title_fontsize,'fontweight','normal', 'Interpreter','tex')
                                            
                                            % format box appearance:
                                            set(gca,'TickDir','out');
                                            box off
                                            ax = axis;
                                            plot(ax(2)*[1,1],ax(3:4),'k','linewidth',1)
                                            plot(ax(1:2),ax(4)*[1,1],'k','linewidth',1)
                                            plot([1,1],ax(3:4),'k','linewidth',1)
                                            plot(ax(1:2),ax(3)*[1,1],'k','linewidth',1)
                                            hold off
                                            
                                            % plot peak threshold line:
                                            if datainfo.peaks_threshold > plotinfo.min_nov_filt(results_count+1)
                                                axes('position', get(gca, 'position'));
                                                pkthr_line = ones(1,mocap_struct_resampled.nFrames) * datainfo.peaks_threshold;
                                                plot(pkthr_line,'--','color',seg_threshold_colour,'LineWidth',1);
                                                set(gca,...
                                                    'YAxisLocation', 'right','color','none',...
                                                    'xlim',[1,mocap_struct_resampled.nFrames],...
                                                    'xtick',[],...
                                                    'ylim',[ plotinfo.min_nov_filt(results_count+1) , plotinfo.max_nov_filt(results_count+1) ],...
                                                    'ytick',datainfo.peaks_threshold,...
                                                    'yticklabel',sprintf('%g\n',datainfo.peaks_threshold),...
                                                    'fontsize',oneplot_ticklabels_fontsize,...
                                                    'TickDir','out')
                                                box off
                                            end
                                            
                                            % display segment numbers:
                                            if segnum_sw
                                                segnum_ind = [1,seg_bounds,mocap_struct_resampled.nFrames];
                                                segnum_y = plotinfo.min_nov_filt(results_count+1) + (plotinfo.max_nov_filt(results_count+1) - plotinfo.min_nov_filt(results_count+1)) * 0.8;
                                                for segnum_count = 1:(length(seg_bounds)+1)
                                                    segnum_thisx = segnum_ind(segnum_count) + (segnum_ind(segnum_count+1) - segnum_ind(segnum_count))/2;
                                                    text(segnum_thisx,segnum_y,num2str(segnum_count),'fontsize',datainfo.segnum_fontsize,...
                                                        'FontWeight','normal','Color',seg_labels_colour,'HorizontalAlignment','Center');
                                                end
                                            end
                                            
                                            drawnow
                                            
                                            % save single-sequence segmentation plot:
                                            if isempty(segplot_save_format) == 0
                                                if iscell(segplot_save_format) == 0
                                                    segplot_save_format = {segplot_save_format};
                                                end
                                                if segnum_sw
                                                    segnum_label = '_numbered';
                                                else
                                                    segnum_label = '';
                                                end
                                                fig_filename = sprintf('%s_%gfps_%s%s_boundaries_%s',...
                                                    file_name,resample_fps,results_label_seg,segnum_label,method_label);
                                                cd(results_path)
                                                for sp_format = segplot_save_format
                                                    if strcmp(sp_format{:},'fig')
                                                        savefig(figh_segplot,[fig_filename,'.fig'])
                                                    elseif strcmp(sp_format{:},'tiff')
                                                        formattype = '-dtiff';
                                                        resolution = '-r300';
                                                        figh_segplot.PaperPositionMode = 'auto';
                                                        print(figh_segplot,[fig_filename,'.tif'],formattype,resolution)
                                                    else
                                                        error('File format for segmentation plot should be either ''fig'' or ''tiff''.')
                                                    end
                                                end
                                            end
                                        end
                                    end
                                end
                                
                                % ..........................................................................
                                % write results record:
                                results_count = results_count + 1;
                                
                                results.distance_measure{results_count,1} = distance_label;
                                results.method{results_count,1} = method_label;
                                results.numvars_values(results_count,:) = [sel_datainfo,mocap_struct_resampled.nMarkers, resample_fps, ...
                                    derivative_order, datainfo.novelty_width_f,datainfo.filter_width_f,...
                                    datainfo.peaks_threshold,datainfo.novelty_width_f/datainfo.filter_width_f];
                               
                                if isempty(economy_sw)==0
                                    results.recording_filename{results_count} = datainfo.file_name;
                                    results.nov_filt{results_count} = nov_filt;
                                    results.seg_bounds{results_count} = seg_bounds;
                                    results.times_values(results_count,1:2) = [ mocap_struct_resampled.nFrames / mocap_struct_resampled.freq , process_stop_time];
                                    results.times_values(results_count,3) = results.times_values(results_count,1)/results.times_values(results_count,2);
                                    % ......................................................................
                                    
                                    % ----------------------------------------------------------------------
                                    % MAKE A STILL OR MOTION PICTURE
                                    
                                    % parameters for preview frame and picture to save:
                                    colour_sw         = 1;          % <--- 0 = black background, 1 = white background
                                    animpar.scrsize   = [400, 400]; % <--- screen size (width, height)
                                    animpar.showmnum  = 0;          % <--- show markers' numbers (0 = no, 1 = yes)
                                    animpar.msize     = 4;          % <--- markers' size
                                    animpar.trl       = 0.3;        % <--- trace length (in seconds, only for preview)
                                    animpar.twidth    = 1;          % <--- trace width
                                    
                                    if colour_sw
                                        animpar.colors = 'wkkkk';
                                    end
                                    if isfield(datainfo,'az_el') == 0 || isempty(datainfo.az_el)
                                        datainfo.az_el = [0, 0];
                                    end
                                    if isfield(datainfo,'limits') && isempty(datainfo.limits) == 0
                                        animpar.limits = datainfo.limits;
                                    else
                                        animpar.limits = [];
                                    end
                                    
                                    % rotate:
                                    if exist('rot_change','var') == 0 ...
                                            || strcmp(rot_change.file_name,datainfo.file_name) == 0 ...
                                            || isequal(rot_change.trim, datainfo.trim) == 0 ...
                                            || isequal(rot_change.select_markers,datainfo.select_markers) == 0 ...
                                            || isequal(rot_change.az_el,datainfo.az_el) == 0
                                        
                                        mocap_struct_rot = mcrotate(mocap_struct_selected,datainfo.az_el(1),[0 0 1]);
                                        mocap_struct_rot = mcrotate(mocap_struct_rot,datainfo.az_el(2),[1 0 0]);
                                        
                                        rot_change.file_name = datainfo.file_name;
                                        rot_change.trim = datainfo.trim;
                                        rot_change.select_markers = datainfo.select_markers;
                                        rot_change.az_el = datainfo.az_el;
                                    end
                                    
                                    % plot preview:
                                    if isempty(preview_frame_time) == 0
                                        animpar.animate  = 0;
                                        mcplotframe(mocap_struct_rot,floor(preview_frame_time*mocap_struct_rot.freq)+1,animpar);
                                    end
                                    
                                    % ..................................................
                                    
                                    if isempty(pic_save_format) == 0
                                        
                                        % parameters for still or motion picture:
                                        preroll           = 0;      % <--- pre-roll time (seconds)
                                        postroll          = 0;      % <--- post-roll time (seconds)
                                        animpar.fps       = 30;     % <--- video frame rate (e.g., 10 for preview; 30 for high quality; 60 for slow motion), also sets resolution for stills
                                        animpar.getparams = 0;      % <--- get parameters to structure 'p' without making picture, especially useful to get optimal limits (0 = no, 1 = yes, pic_save_format should not be empty)
                                        
                                        segments_n = length(seg_bounds)+1;
                                        
                                        if isempty(all_selected_segments)
                                            selected_segments = [];
                                        else
                                            if isnumeric(all_selected_segments)
                                                if  max(all_selected_segments) > segments_n
                                                    error('Selected segments to make picture exceed computed segments = %i',segments_n)
                                                else
                                                    selected_segments = all_selected_segments;
                                                end
                                            else
                                                selected_segments = [1:segments_n];
                                            end
                                        end
                                        if iscell(pic_save_format) == 0
                                            pic_save_format = {pic_save_format};
                                        end
                                        
                                        % make stills or motion picture:
                                        seg_bounds_plus = [1, floor(seg_bounds * (mocap_struct_rot.freq / resample_fps)),  mocap_struct_rot.nFrames];
                                        if isempty(selected_segments)
                                            vloop_laps = 1;
                                        else
                                            vloop_laps = length(selected_segments);
                                        end
                                        for pic_format = pic_save_format
                                            if strcmp(pic_format{:},'tiff') || strcmp(pic_format{:},'fig')
                                                animpar.animate     = 0;      % 0 = last frame of segment, 1 = video
                                                animpar.trl         = 'full'; % trace length (in seconds, 'full' = trace whole segment)
                                                animpar.frameformat = pic_format{:};
                                            elseif strcmp(pic_format{:},'mpeg4')
                                                animpar.animate     = 1;
                                                animpar.trl         = video_trace_length;
                                                animpar.videoformat = pic_format{:};
                                            end
                                            animpar.createframes = abs(1-animpar.animate);
                                            if animpar.trl == 0
                                                trace_label = 'notrace';
                                            elseif isnumeric(animpar.trl) == 0 && strcmp(animpar.trl,'full')
                                                trace_label = 'fulltrace';
                                            else
                                                if animpar.animate
                                                    trace_label = sprintf('trace%g',animpar.trl);
                                                end
                                            end
                                            results_label_az_el = sprintf('az%g_el%g_%s',datainfo.az_el(1),datainfo.az_el(2),trace_label);
                                            pic_path = sprintf('%s_%gfps_%s_%s_%s_%s',file_name,resample_fps,results_label_seg,results_label_az_el,upper(pic_format{:}),method_label);
                                            cd(results_path);
                                            if exist(pic_path,'dir') ~= 7 && isempty(selected_segments) == 0 && animpar.getparams == 0
                                                mkdir(pic_path);
                                            end
                                            preroll_f = preroll * mocap_struct_rot.freq;
                                            postroll_f = postroll * mocap_struct_rot.freq;
                                            for picseg_count = 1:vloop_laps
                                                if isempty(selected_segments) == 0 && sum(selected_segments==0)==0 % generate pictures for specified segments
                                                    mocap_struct_segment = ...
                                                        mctrim(mocap_struct_rot,seg_bounds_plus(selected_segments(picseg_count))-preroll_f,seg_bounds_plus(selected_segments(picseg_count)+1)+postroll_f,'frame');
                                                    animpar.output = sprintf('%s/%s_%gfps_%s_s%d_%s',pic_path,file_name,resample_fps,results_label_seg,selected_segments(picseg_count),results_label_az_el);
                                                else % generate picture for full data length
                                                    mocap_struct_segment = mocap_struct_rot;
                                                    animpar.output = sprintf('%s_full_az%g_el%g_%s',file_name,datainfo.az_el(1),datainfo.az_el(2),trace_label);
                                                end
                                                if strcmp(animpar.trl,'full') &&  animpar.animate == 0
                                                    animpar.frames = mocap_struct_segment.nFrames;
                                                else
                                                    animpar.frames = [];
                                                end
                                                p  = mcanimate(mocap_struct_segment,animpar);
                                            end
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end

% ----------------------------------------------------------------------
% PLOT AND SAVE ALL RESULTS

if plotsaveallresults_sw
    
    % save results:
    cd(results_path)
    pause(1) % ensures that there are no two equal time stamps, since it is rounded to seconds
    results_filename = sprintf('results_%d_%d_%d_%d_%d_%d',fix(clock)); % results_year_month_day_hour_minute_second
    save(results_filename,'results');

    % plot all sequences:
    seg_threshold_width = 2;      % <--- segmentation threshold line width
    alph_label_vo       = -0.023; % <--- vertical offset for alphabetic labels
    alph_label_ho       = 0.09;   % <--- horizontal offset for alphabetic labels
    alph_label_fs       = 18;     % <--- font size for alphabetic labels
    
    allsegfig_monpos = selected_monpos;
    allsegfig_monpos(1) = selected_monpos(3)/4;
    allsegfig_monpos(3) = selected_monpos(3)/2;

    figh_allsegplot = figure('Position',allsegfig_monpos);
    n_results = length(results.nov_filt);
    trial_index = 1;
    
    for subplot_index = 1:n_results
        
        sp_ax{subplot_index} =  subplot(n_results,1,subplot_index);
        
        % plot filtered novelty score:
        plot(results.nov_filt{subplot_index},'color',seg_curve_colour)
        set(gca,...
            'xlim',[1,length(results.nov_filt{subplot_index})],...
            'xtick',[],...
            'xticklabel',[],...
            'ylim',[ plotinfo.min_nov_filt(subplot_index) , plotinfo.max_nov_filt(subplot_index) ],...
            'ytick',[],...
            'yticklabel',[])
        hold on
        
        % plot computed boundaries:
        rep_peaks = repmat(results.seg_bounds{subplot_index},2,1);
        lines = repmat(get(gca,'ylim'),size(rep_peaks,2),1)';
        line(rep_peaks,lines,'Color',seg_boundaries_colour,'linewidth',3,'linestyle','-');
        
        % format box appearance:
        set(gca,'TickDir','out');
        box off
        ax = axis;
        plot(ax(2)*[1,1],ax(3:4),'k','linewidth',1)
        plot(ax(1:2),ax(4)*[1,1],'k','linewidth',1)
        plot([1,1],ax(3:4),'k','linewidth',1)
        plot(ax(1:2),ax(3)*[1,1],'k','linewidth',1)
        hold off
                
        % plot peak threshold line:
        if results.numvars_values(subplot_index,7) > 0
            pkthr_line = ones(1,length(results.nov_filt{subplot_index})) * results.numvars_values(subplot_index,7);
            hold on
            plot(pkthr_line,'--','color',seg_threshold_colour,'LineWidth',seg_threshold_width); % threshold line
            hold off
            threshtickyaxis_pos = get(gca, 'position');
            threshtickyaxis_pos(1) = threshtickyaxis_pos(1) + threshtickyaxis_pos(3);
            threshtickyaxis_pos(3) = threshtickyaxis_pos(3) * 0.01;
            axes('position',threshtickyaxis_pos);
            yaxis_line = ones(1,2) * results.numvars_values(subplot_index,7);
            plot(yaxis_line,'-','color',seg_threshold_colour,'LineWidth',seg_threshold_width); % tick
            ytick_ylim = [ plotinfo.min_nov_filt(subplot_index) , plotinfo.max_nov_filt(subplot_index) ];
            set(gca,...
                'Color','none','XColor','none','YColor','none',...
                'ylim',[1,2],...
                'ylim',ytick_ylim)
            text(3,results.numvars_values(subplot_index,7),sprintf('\\theta'),'fontsize',alph_label_fs,'fontweight','bold')
            box off
        end
        
        % alphanumerical labels to the left of the figure, in front of each panel
        alphabet = 'abcdefghijklmnopqrstuvwxyz';
        if subplot_index > 1
            if results.numvars_values(subplot_index,1) == results.numvars_values(subplot_index-1)
                trial_index = trial_index + 1;
            else
                trial_index = 1;
            end
        end
        results.numvars_values(subplot_index,1);
        row_labels_axes{subplot_index} = axes('Units','normalized','Position',...
            [ alph_label_ho, sp_ax{subplot_index}.Position(2)+alph_label_vo, 0, sp_ax{subplot_index}.Position(4)/2 ],'Visible','off');
        set(get(row_labels_axes{subplot_index},'Title'),'Visible','on')
        row_labels{subplot_index} = title(sprintf('%d%s',results.numvars_values(subplot_index,1),alphabet(trial_index)),'fontsize',alph_label_fs,...
            'VerticalAlignment','bottom','HorizontalAlignment','center','Interpreter','none','FontWeight','Normal');
    end
    
    % save multi-sequence segmentation plot:
    if isempty(segplot_save_format) == 0
        if iscell(segplot_save_format) == 0
            segplot_save_format = {segplot_save_format};
        end
        cd(results_path)
        for sp_format = segplot_save_format
            if strcmp(sp_format{:},'fig')
                savefig(figh_allsegplot,[results_filename,'.fig'])
            elseif strcmp(sp_format{:},'tiff')
                formattype = '-dtiff';
                resolution = '-r300';
                figh_allsegplot.PaperPositionMode = 'auto';
                print(figh_allsegplot,[results_filename,'.tif'],formattype,resolution)
            else
                error('File format for segmentation plot should be either ''fig'' or ''tiff''.')
            end
        end
    end
end

program_stop_time = toc(program_start_time);
timer_hours = floor(program_stop_time / 3600);
remainder_seconds = program_stop_time - timer_hours * 3600;
timer_minutes = floor(remainder_seconds / 60);
timer_seconds = remainder_seconds - timer_minutes * 60;
disp(sprintf('program time: %2d hours, %2d minutes, %f seconds\n', timer_hours, timer_minutes, timer_seconds))
