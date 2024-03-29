%% Script developed to read and modify Abaqus input files %%
% E. Martínez-Pañeda and M. Muñiz-Calvente
% mail@empaneda.com and miguelmunizcalvente@gmail.com

% If using this script for research or industrial purposes, please cite:
% G. Papazafeiropoulos, M. Muñiz-Calvente, E. Martínez-Pañeda. 
% Abaqus2Matlab: a suitable tool for finite element post-processing.
% Advances in Engineering Software 105, 9-16 (2017)
% doi:10.1016/j.advengsoft.2017.01.006

% CPE8 ELEMENT

clear variables
clc

% Open the file
fid = fopen('Job-1.inp','r');
%Pasar el archivo a caracteres, detecta cambio de línea como delimitador checkbox_y
%si encuentra una línea en blanco la sustituye por comillas
file=textscan(fid,'%s','Delimiter','\n');
%fclose(fid);
%Para que quede en una columna del tipo Cell Array
flines=(file{1});
%Donde acabamos de leer
lif=length(flines); %última linea del fichero
%fclose(fid);

% Search the lines containing *ELEMENT
FFLINES=upper(flines);
Lin_element_cells=strfind(FFLINES,['*' upper('Element')]); 
Lin_element_zeros=cellfun(@isempty,Lin_element_cells); 
start_elements=find(Lin_element_zeros==0); %numbers of the lines that contains *Element

Lin_order_cells=strfind(FFLINES,'*');
Lin_order_zeros=cellfun(@isempty,Lin_order_cells); 
Lin_orders=find(Lin_order_zeros==0); %numbers of the lines that contains *Element

for i=1:length(start_elements)
 end_elements(i)=Lin_orders(find(Lin_orders==start_elements(i))+1);
 p=1;
 str_a=FFLINES{start_elements(i)}(~isspace(FFLINES{start_elements(i)}));%eliminate spaces
 flag=0;
 if length(str_a)==(length('Element')+1)
  flag=1;    
 else
  if str_a(length('Element')+2)==','
   flag=1;
  end
 end

 if flag==1
  for b=(start_elements(i)+1):(end_elements(i)-1) 
   Element{i}(p,:)=strread(FFLINES{b}, '%f', 'delimiter', ',');
   p=p+1;
  end
 end
end

LE=Element{1}; % Element numbering and connectivity
[NUME, ~] = size(LE);
Order=floor(log10(NUME));
offset=10*10^Order;
LEnew=zeros(NUME,9);
for i=1:NUME
 LEnew(i,:)= [(offset+i) LE(i,2) LE(i,3) LE(i,4) LE(i,5) LE(i,6) LE(i,7) LE(i,8) LE(i,9)];   
end    
 for i=1:length(Element)
  % A new file is created so as to avoid overwriting the previous one
  NameFile=strcat('Job-1a.inp');
  % We open the new input file and start making changes
  dlmwrite(NameFile,LEnew,'precision', 8);
 end

fclose all;