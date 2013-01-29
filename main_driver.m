%% The driver program to run all the component in this folder
% Written by Ge Jin,jinwar@gmail.com

% automatic select the window range
run_autowinpick

% making cross-correlation measurement
gsdfmain

% calculating eikonal tomography for each event
eikonal_eq

% stacking the result
stack_phv
