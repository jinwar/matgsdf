%% The driver program to run all the component in this folder
% Written by Ge Jin,jinwar@gmail.com

% if using sac dataset
%sac2eventmat

% download the data using 
data_download

% clean up multiple or close events
cleanup_events

% automatic select the window range
run_autowinpick

% making cross-correlation measurement
gsdfmain

% calculating eikonal tomography for each event
eikonal_eq

% stacking the result
stack_phv

% apply amplitude correction
helmholtz_eq

% stack the result of helmholtz
stack_helm

% export the result into xyz format
make_xyz
