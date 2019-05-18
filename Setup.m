% Run all settings and save figures
clear all 
clc
%cdgea;
addpath('tools')
%addpath(genpath('F:\jamesjcai-scGEApp-65b47a5')) 
%% import data
%import_gene_lists

%import_Liver_data;

%% generate colors and symbols
cols_symbs

%% classify cell types and generate figures
celltype_classification('8k','merged',1);

close all;

