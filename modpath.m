function modpath()

this_file_name=mfilename('fullpath');
this_dir_name=fileparts(this_file_name);
addpath(fullfile(this_dir_name,'toolbox'));
addpath(genpath(fullfile(this_dir_name,'toolbox','tmt_rel_114')));
addpath(genpath(fullfile(this_dir_name,'toolbox','snippeter')));
%addpath(genpath(fullfile(this_dir_name,'taylor_matlab_toolbox/113')));
%addpath(genpath(fullfile(this_dir_name,'groundswell/groundswell_rel_118')));
parent_dir_name=fileparts(this_dir_name);
addpath(degit(genpath(fullfile(parent_dir_name,'ax_repo'))));
addpath(degit(genpath('~/groundswell/repo')));

end
