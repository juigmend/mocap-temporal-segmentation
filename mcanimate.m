function p = mcanimate(d, p, proj)
% Creates animation of mocap data and saves it to file (.avi or .mpeg-4) or as consecutive
% frames (.png). Matlab's VideoWriter function is used to create the video
% file.
% COMPATIBILITY NOTES (v. 1.5): The 'folder'-field (animpar structure, v. 1.4) has been changed to 
% 'output' and is used as file name for the animation (and stored to the current directory) 
% or as folder name in case frames are to be plotted. 
% Please use the function without the projection input argument, but specify 
% it in the animation structure instead.
%
% syntax
% p = mcanimate(d);
% p = mcanimate(d, p);
% 
% input parameters
% d: MoCap data structure
% p: animpar structure (optional)
% [depricated: proj: projection used: 0 = orthographic (default), 1 = perspective
%                    this flag is supposed to be set in the animation parameter stucture instead]
% p.frames : frames to plot (optional), if only one frame, then a folder will not be created  % (JIMG 2018/10/25)
%
% output
% p: animpar structure used for plotting the frames
% 
% examples
% mcanimate(d, par);
%
% comments
% If the animpar structure is not given as input argument, the function
% creates it by calling the function mcinitanimpar and setting the .limits field 
% of the animpar structure automatically so that all the markers fit into all frames.
% If the par.pers field (perspective projection) is not given, it is created internally for backwards
% compatibility. For explanation of the par.pers field, see help mcinitanimpar                                                                                
%
% see also
% mcplotframe, mcinitanimpar
%
% Part of the Motion Capture Toolbox, Copyright 2008, 
% University of Jyvaskyla, Finland

% Lines changed by Juan Ignacio Mendoza Garay marked with (JIMG year/month/day) where year/month/day is the date of the change.

if nargin>2
    disp([10, 'Please note that, from MCT version 1.5, the perspective projection flag is set in the field "perspective" in the animation parameters. Please adapt your code accordingly.', 10])
    p.perspective=proj;
end

if nargin<2
    p = mcinitanimpar;
end

% resample for p.fps fps movie
if size(d.data,1) > 1                                                           % (JIMG 2019/1/30)
    d2 = mcresample(d, p.fps);
else                                                                            % (JIMG 2019/1/30)
    d2 = d;                                                                     % (JIMG 2019/1/30)
end                                                                             % (JIMG 2019/1/30)

p.animate = 1;

if isfield(p,'frames') == 0 || isempty(p.frames) || p.frames == 0                % (JIMG 2018/10/25)
    frames = 1:d2.nFrames;                                                       % (JIMG 2018/10/25)
else                                                                             % (JIMG 2018/10/25)
    frames = floor(p.fps * p.frames/d.freq);                                     % (JIMG 2019/1/29)
end                                                                              % (JIMG 2018/10/25)

if frames == 0                                                                   % (JIMG 2019/1/30)
    frames = 1;                                                                  % (JIMG 2019/1/30)
end                                                                              % (JIMG 2019/1/30)

% p = mcplotframe(d2,1:d2.nFrames,p); % output parameter added 240608  % commented (JIMG 2018/10/25)
p = mcplotframe(d2,frames,p);                                                    % (JIMG 2018/10/25)


