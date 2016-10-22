
% This is one way of doing indirect inference:
% open up stata, run the pvar command, then import the coefs

% problem: it takes one second per iteration. perhaps too slow
% in comparison, the non-indirect-inference approach takes 0.05 seconds per
% iteration

tic
%% Run stata
system('"C:\Program Files (x86)\Stata13\StataMP-64.exe" /e do C:\Users\pedm\Documents\Research\ProductivityProject\Do\pvar_matlab_data.do');
toc

%% Read in the pvar coefs
input = caseread('pvar_coefs_matlab.txt');
split = strsplit(input(2,:), '	');
pvar_coefs = str2double(split(2:end-1))

%% This should always be 9
length(pvar_coefs)
