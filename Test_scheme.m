clear all
close all
clc

addpath('Resources') 
%% Generate new folder
folder_name='Results/Prova'; %Name for the results folder: it should be named after the kind of test performed

currDate = datestr(datetime,30);
mkdir(folder_name,currDate);
results_folder=strcat(folder_name,'/',currDate);
%All the figures we want to save should be added to the array of figures
%through figs(1)=figure(1); figs(2)=figure(2);

%% Do whatever you need to do
a=10;
b=linspace(1,10,1000);
c=zeros(100,100);

figs(1)=figure(1);
plot(1:100,1:100);
hold on
figs(2)=figure(2);
plot(1:10,cos(1:10));
figs(1)=figure(1);
plot(1:100,cos(1:100));

min_dist_roots(b)

%% Save workspace and figures to the folder
savefig(figs,strcat(results_folder,'/figures.fig'),'compact');
clear figs
save(strcat(results_folder,'/workspace'));