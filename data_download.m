clear;

javaaddpath('IRIS-WS-2.0.6.jar');

setup_parameters;

if exist('timestamp.mat','file')
	temp = load('timestamp.mat');
	start_time = temp.start_time;
else
	start_time = parameters.start_time;
end


