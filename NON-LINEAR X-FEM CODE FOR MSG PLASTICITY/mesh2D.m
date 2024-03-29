% Function to read nodal coordinates and element connectivity from text
% files. A commercial FE package or any other mesh generator can be used to
% create 2 files: nodes.txt (nodal coordinates) and conec.txt (element
% connectivity). Current version: Q4 elements

% If using this code for research or industrial purposes please cite:
% E. Martínez-Pañeda, S. Natarajan, S. Bordas. 
% Gradient plasticity crack tip characterization by means of the extended
% finite element method. Computational Mechanics (2017)
% doi:10.1007/s00466-017-1375-6

function [XY,LE]=mesh2D

fileID = fopen('nodes.txt','r');
formatSpec = '%d %f %f';
sizeA = [3 Inf];
Ab = fscanf(fileID,formatSpec,sizeA);
fclose(fileID);
XY=Ab';
XY(:,1)=[];

fileID = fopen('conec.txt','r');
formatSpec = '%d %d %d %d %d';
sizeA = [5 Inf];
Ac = fscanf(fileID,formatSpec,sizeA);
LE=Ac';
LE(:,1)=[];

% Emilio Martínez-Pañeda (mail@empaneda.com)