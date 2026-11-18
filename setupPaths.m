function root = setupPaths()
%SETUPPATHS  Put the blueprint package and the host simulator's Patient
%   class on the MATLAB path. Call once per session.
root = fileparts(mfilename('fullpath'));
addpath(fullfile(root,'models'), fullfile(root,'core'), ...
        fullfile(root,'experiments'));
hostInit = fullfile(fileparts(root),'initialization');
if exist(hostInit,'dir'), addpath(hostInit); end
% fall back to the bundled copy of Patient when the package is distributed
% on its own, without the host simulator around it
if isempty(which('Patient')), addpath(fullfile(root,'compat')); end
end
