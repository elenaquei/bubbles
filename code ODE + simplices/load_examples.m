% global variables used by the library - reset to specified interest
global debug
debug = 0;

global nu
nu = 1.1;

global use_intlab
use_intlab = 0;

global talkative
talkative = 1;

global RAD_MAX
RAD_MAX = 10^-4;

global norm_weight
norm_weight = [];

global Display
Display = 0;

% start Intlab
startintlab

% find location of library
path = mfilename('fullpath');
[filepath,name,ext] = fileparts(path);

% find if if uses \ or / 
if any(filepath == '\')
    bar = '\';
else
    bar = '/';
end

% add to path the interesting folders
addpath(append(filepath, bar, 'branch Hopf example'))
addpath(append(filepath, bar, 'branch ODE templates'))
addpath(append(filepath, bar, 'simplex DDE Hopf examples'))
addpath(append(filepath, bar, 'simplex ODE Hopf examples'))
