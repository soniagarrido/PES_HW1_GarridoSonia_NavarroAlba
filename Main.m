% This program solves a convection-diffusion problem 
% 

clear, close all, home

tic

global diffusion  h 

disp(' ')
disp('No source term is considered');

diffusion = input('Diffusion coefficient = ');
disp(' ')

disp(' Select the problem type');
disp('1: 2D_quad_linear');
disp('2: 2D_quad_quad');
disp('3: 2D_tri_linear');
disp('4: 2D_tri_quad');
disp('5: 3D_quad_linear');
disp('6: 3D_quad_quad');
disp('7: 3D_tri_linear');
disp('8: 3D_tri_quad');
Problem=input('Problem type   ');
if Problem ~= [1 2 3 4 5 6 7 8]
    disp ('ERROR: Select an existing problem type');
    clear all
    return
end

if Problem==1
X1=load('nodes_2D_quad_linear');
X=X1(:,2:3);
T1=load('elm_2D_quad_linear');
T=T1(:,2:5);

ngaus = 4;
% Quadrature,Shape Functions
[n,wpg,pospg,N,dNdxi] = C2D4 ;

% SYSTEM RESULTING OF DISCRETIZING THE WEAK FORM
[K,f] = CreateMatrix(X,T,pospg,wpg,N,dNdxi);

% BOUNDARY CONDITIONS 

    % nodes on which solution is u=1
    nodesDir1 = load ('inlet.dat') ;
    % nodes on which solution is u=0
    nodesDir0 = load('outlet.dat');
    % Boundary condition matrix
    C = [nodesDir1, ones(length(nodesDir1),1);
         nodesDir0, zeros(length(nodesDir0),1)];
elseif Problem==2
    X1=load('nodes_2D_quad_quad');
    X=X1(:,2:3);
    T1=load('elm_2D_quad_quad');
    T=T1(:,2:9);
    % Number of gauss points
    ngaus = 9;
    % Quadrature,Shape Functions
    [n,wpg,pospg,N,dNdxi] = C2D8 ;
    pospg=pospg';
    % SYSTEM RESULTING OF DISCRETIZING THE WEAK FORM
    [K,f] = CreateMatrix(X,T,pospg,wpg,N,dNdxi);
    % BOUNDARY CONDITIONS 

    % nodes on which solution is u=1
    nodesDir1 = [1,  14,  75,  76,  77,  78,  79,  80, 480, 501, 505, 607, 641, 647, 714]';
    % nodes on which solution is u=0
    nodesDir0 = [4,   5,   6,   7,   8,   9,  10,  11,  12, 283, 291, 487, 488, 489, 496, 728, 731]';
    % Boundary condition matrix
    C = [nodesDir1, ones(length(nodesDir1),1);
         nodesDir0, zeros(length(nodesDir0),1)];
elseif Problem==3
    X1=load('nodes_2D_tri_linear');
    X=X1(:,2:3);
    T1=load('elem_2D_tri_linear');
    T=T1(:,2:4);
    % NUMERICAL INTEGRATION
    % Number of gauss points
    ngaus = 1;
    % Quadrature,Shape Functions
    [n,wpg,pospg,N,dNdxi] = C2D3 ;
    N=N';
    dNdxi=dNdxi';
    % SYSTEM RESULTING OF DISCRETIZING THE WEAK FORM
    [K,f] = CreateMatrix(X,T,pospg,wpg,N,dNdxi);
    % BOUNDARY CONDITIONS 
    % nodes on which solution is u=1
    nodesDir1 = load ('inlet.dat') ;
    % nodes on which solution is u=0
    nodesDir0 = load('outlet.dat');
    % Boundary condition matrix
    C = [nodesDir1, ones(length(nodesDir1),1);
         nodesDir0, zeros(length(nodesDir0),1)];

elseif Problem==4
    X1=load('nodes_2D_tri_quad');
    X=X1(:,2:3);
    T1=load('elem_2D_tri_quad');
    T=T1(:,2:7);
    % NUMERICAL INTEGRATION
    % Number of gauss points
    ngaus = 3;
    % Quadrature,Shape Functions
    [n,wpg,pospg,N,dNdxi] = C2D6 ;
    pospg=pospg';


    % SYSTEM RESULTING OF DISCRETIZING THE WEAK FORM
    [K,f] = CreateMatrix(X,T,pospg,wpg,N,dNdxi);

    % BOUNDARY CONDITIONS 
    % nodes on which solution is u=1
    nodesDir1 = [1,  14,  75,  76,  77,  78,  79,  80, 453, 462, 696, 698, 701, 705, 843]';
    % nodes on which solution is u=0
    nodesDir0 = [4,   5,   6,   7,   8,   9,  10,  11,  12, 263, 465, 466, 468, 470, 474, 478,711]';
    % Boundary condition matrix
    C = [nodesDir1, ones(length(nodesDir1),1);
         nodesDir0, zeros(length(nodesDir0),1)];
elseif Problem==5
    X1=load('nodes_3D_quad_lin.txt');
    X=X1(:,2:4);
    T1=load('elem_3D_quad_lin.txt');
    T=T1(:,2:9);
    % NUMERICAL INTEGRATION
% Number of gauss points
ngaus = 8;
% Quadrature,Shape Functions
[n,wpg,pospg,N,dNdxi] = C3D8;
pospg=pospg';

% SYSTEM RESULTING OF DISCRETIZING THE WEAK FORM
[K,f] = CreateMatrix_3D(X,T,pospg,wpg,N,dNdxi,ngaus);

% BOUNDARY CONDITIONS 
    % nodes on which solution is u=1
    nodesDir1 = [1,   3,  24,  25,  26,  27,  28,  29, 260, 262, 283, 284, 285, 286, 287, 288, 519, 521, 542, 543, 544, 545, 546, 547]';
    % nodes on which solution is u=0
    nodesDir0 = [5,   6,   7,   8,   9,  10,  11,  12,  13, 264, 265, 266, 267, 268, 269, 270, 271, 272, 523, 524, 525, 526, 527, 528, 529, 530, 531]';
    % Boundary condition matrix
    C = [nodesDir1, ones(length(nodesDir1),1);
         nodesDir0, zeros(length(nodesDir0),1)];
elseif Problem==6
    X1=load('nodes_3D_quad_quad.txt');
    X=X1(:,2:4);
    T1=load('elem_3D_quad_quad.txt');
    T=T1(:,2:21);
    % NUMERICAL INTEGRATION
% Number of gauss points
ngaus = 27;
% Quadrature,Shape Functions
[n,wpg,pospg,N,dNdxi] = C3D20;
pospg=pospg';
% N=N';
% dNdxi=dNdxi';

% SYSTEM RESULTING OF DISCRETIZING THE WEAK FORM
[K,f] = CreateMatrix_3D(X,T,pospg,wpg,N,dNdxi,ngaus);

% BOUNDARY CONDITIONS 
    % nodes on which solution is u=1
    nodesDir1 = [1,    3,   24,   25,   26,   27,   28,   29,  260,  262,  283,  284,  285,  286,  287,  288, 519,  521,  542,  543,  544,  545,  546,  547,  903,  906,  909,  911,  975,  982,  984,  986, 992,  997,  998,  999, 1000, 1001, 1449, 1456, 1457, 1459, 1464, 1466, 1772, 1773, 2070, 2073, 2075, 2119, 2121, 2123, 2129, 2130, 2131, 2132, 2414, 2418, 2422, 2424, 2606]';
    % nodes on which solution is u=0
    nodesDir0 = [5,    6,    7,    8,    9,   10,   11,   12,   13,  264,  265,  266,  267,  268,  269,  270, 271,  272,  523,  524,  525,  526,  527,  528,  529,  530,  531,  839,  844,  848,  849,  852, 855,  858,  860,  862,  865, 1413, 1420, 1422, 1423, 1424, 1431, 1432, 1434, 1974, 1977, 1978, 1979, 1982, 1983, 1984, 2031, 2035, 2036, 2038, 2041, 2043, 2046, 2395, 2397, 2398, 2402, 2403, 2405, 2717, 2718, 2720, 2721]';
    % Boundary condition matrix
    C = [nodesDir1, ones(length(nodesDir1),1);
         nodesDir0, zeros(length(nodesDir0),1)];
elseif Problem==7
    X1=load('nodes_3D_tri_lin.inp');
    X=X1(:,2:4);
    T1=load('elem_3D_tri_lin.txt');
    T=T1(:,2:5);
    % NUMERICAL INTEGRATION
% Number of gauss points
ngaus = 1;
% Quadrature,Shape Functions
[n,wpg,pospg,N,dNdxi] = C3D4;
pospg=pospg';
N=N';
dNdxi=dNdxi';

% SYSTEM RESULTING OF DISCRETIZING THE WEAK FORM
[K,f] = CreateMatrix_3D(X,T,pospg,wpg,N,dNdxi,ngaus);

% BOUNDARY CONDITIONS 
    % nodes on which solution is u=1
    nodesDir1 = [5,   6,  11,  12,  54,  98,  99, 100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 207, 208, 209, 210, 211, 212]' ;
    % nodes on which solution is u=0
    nodesDir0 = [7,   8,   9,  10,  64,  65,  66,  67,  68,  69,  70,  71,  72,  73,  74,  75,  76,  77,  78,  79, 191, 192, 193, 194, 195, 196, 197]';
    % Boundary condition matrix
    C = [nodesDir1, ones(length(nodesDir1),1);
         nodesDir0, zeros(length(nodesDir0),1)];
elseif Problem==8
    X1=load('nodes_3D_tri_quad.txt');
    X=X1(:,2:4);
    T1=load('elem_3D_tri_quad.inp');
    T=T1(:,2:11);
    % NUMERICAL INTEGRATION
% Number of gauss points
ngaus = 4;
% Quadrature,Shape Functions
[n,wpg,pospg,N,dNdxi] = C3D10;
pospg=pospg';
% N=N';
% dNdxi=dNdxi';

% SYSTEM RESULTING OF DISCRETIZING THE WEAK FORM
[K,f] = CreateMatrix_3D(X,T,pospg,wpg,N,dNdxi,ngaus);

% BOUNDARY CONDITIONS 
    % nodes on which solution is u=1
    nodesDir1 = [5,    6,   11,   12,   54,   98,   99,  100,  101,  102,  103,  104,  105,  106,  107,  108, 109,  110,  207,  208,  209,  210,  211,  212, 1007, 1012, 1017, 1088, 1102, 1103, 1113, 1204, 1205, 1316, 1318, 1319, 1352, 1720, 1724, 1725, 1813, 1815, 1817, 1852, 1854, 1978, 2094, 2095, 2133, 2259, 2260, 2269, 2271, 2273, 2385, 2431, 2515, 2516, 2833, 2835, 2975, 3039, 3190, 3462,  3463, 3511, 3512, 3683, 3691, 3877, 4097, 4414, 4428, 4429, 4521]' ;
    % nodes on which solution is u=0
    nodesDir0 = [7,    8,    9,   10,   64,   65,   66,   67,   68,   69,   70,   71,   72,   73,   74,   75, 76,   77,   78,   79,  191,  192,  193,  194,  195,  196,  197, 1033, 1133, 1136, 1141, 1168,1174, 1239, 1241, 1243, 1248, 1311, 1315, 1486, 1497, 1506, 1668, 1669, 1761, 1763, 2810, 2857, 2859, 2860, 3468, 3469, 3697, 3936, 3946, 3948, 4114, 4116, 4117, 4194, 4196, 4228, 4229, 4249, 4452, 4455, 4456, 4461, 4463, 4464, 4465, 4466, 4467, 4470, 4471, 4472, 4474, 4475, 4477, 4484, 4485, 4486, 4489, 4490, 4507]';
    % Boundary condition matrix
    C = [nodesDir1, ones(length(nodesDir1),1);
         nodesDir0, zeros(length(nodesDir0),1)];
end

ndir = size(C,1);
neq  = size(f,1);
A = zeros(ndir,neq);
A(:,C(:,1)) = eye(ndir);
b = C(:,2);


% SOLUTION OF THE LINEAR SYSTEM
% Entire matrix
Ktot = [K A';A zeros(ndir,ndir)];
ftot = [f;b];

sol = Ktot\ftot;
Temp = sol(1:neq);
multip = sol(neq+1:end);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % POSTPROCESS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Getting geometry from abaqus to export it to ensight
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%% LOADING FILES OF NODES AND ELEM %%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%% WRITING HEADING FOR VTK FILE

if Problem==1
    
X1=load('nodes_2D_quad_linear');
X=X1(:,2:3);
T1=load('elm_2D_quad_linear');
T=T1(:,2:5);
nnode=length(X);
nelem=length(T);
for i=1:1:nnode
    X(i,3)=0;
end
	% printing heading to file
f=fopen('MyParaviewFile_2D_quad_linear.vtk','w');
fprintf(f,'# vtk DataFile Version 1.0\n');
fprintf(f,'ECM-CELL DIFFUSION-MECHANICS\n');
fprintf(f,'ASCII\n');
fprintf(f,'\n');
fprintf(f,'DATASET UNSTRUCTURED_GRID\n');
fprintf(f,'%s %8i %s\n','POINTS', nnode ,'float');

%%%%%%%%%%%%%%%%%%%%%% WRITING COORDINATES OF NODES %%%%%%%%%%%%%%%%%%%%%%%%%

	% printing coordinates
for i1=1:nnode
           fprintf(f,'%14.8E %14.8E %14.8E\n',X(i1,1:3));
end
fprintf(f,'\n');

%%%%%%%%%%%%%%%%%%%%%% WRITING CONNECTIVITY OF NODES %%%%%%%%%%%%%%%%%%%%%%%%%

fprintf(f,'%s %8i %8i\n','CELLS', nelem , nelem*5);
for i1=1:nelem
                
            fprintf(f,'%4i  %10i  %10i  %10i %10i\n',4,T(i1,1)-1,T(i1,2)-1,T(i1,3)-1,T(i1,4)-1);
end
fprintf(f,'\n');
fprintf(f,'%s %8i\n','CELL_TYPES', nelem);
for i1=1:nelem
            fprintf(f,' %4i ', 9);
end
fprintf(f,'\n');          

%%%%%%%%%%%%%%%%%%%%%% WRITING VARIABLES %%%%%%%%%%%%%%%%%%%%%%%%%
fprintf(f,'%s %8i\n','POINT_DATA', nnode);
fprintf(f,'SCALARS Ux float\n');
fprintf(f,'LOOKUP_TABLE default\n');
for i1=1:nnode
           fprintf(f,'%14.8E\n',Temp(i1) );
end
elseif Problem==2
    X1=load('nodes_2D_quad_quad');
X=X1(:,2:3);
T1=load('elm_2D_quad_quad');
T=T1(:,2:9);
nnode=length(X);
nelem=length(T);
for i=1:1:nnode
    X(i,3)=0;
end
	% printing heading to file
f=fopen('MyParaviewFile_2D_quad_quad.vtk','w');
fprintf(f,'# vtk DataFile Version 1.0\n');
fprintf(f,'ECM-CELL DIFFUSION-MECHANICS\n');
fprintf(f,'ASCII\n');
fprintf(f,'\n');
fprintf(f,'DATASET UNSTRUCTURED_GRID\n');
fprintf(f,'%s %8i %s\n','POINTS', nnode ,'float');

%%%%%%%%%%%%%%%%%%%%%% WRITING COORDINATES OF NODES %%%%%%%%%%%%%%%%%%%%%%%%%

	% printing coordinates
for i1=1:nnode
           fprintf(f,'%14.8E %14.8E %14.8E\n',X(i1,1:3));
end
fprintf(f,'\n');

%%%%%%%%%%%%%%%%%%%%%% WRITING CONNECTIVITY OF NODES %%%%%%%%%%%%%%%%%%%%%%%%%

fprintf(f,'%s %8i %8i\n','CELLS', nelem , nelem*9);
for i1=1:nelem                
            fprintf(f,'%4i  %10i  %10i  %10i %10i %10i %10i %10i %10i\n',8,T(i1,1)-1,T(i1,2)-1,T(i1,3)-1,T(i1,4)-1,T(i1,5)-1,T(i1,6)-1,T(i1,7)-1,T(i1,8)-1);
end
fprintf(f,'\n');
fprintf(f,'%s %8i\n','CELL_TYPES', nelem);
for i1=1:nelem
            fprintf(f,' %4i ', 23);
end
fprintf(f,'\n');          

%%%%%%%%%%%%%%%%%%%%%% WRITING VARIABLES %%%%%%%%%%%%%%%%%%%%%%%%%
fprintf(f,'%s %8i\n','POINT_DATA', nnode);
fprintf(f,'SCALARS Ux float\n');
fprintf(f,'LOOKUP_TABLE default\n');
for i1=1:nnode
           fprintf(f,'%14.8E\n',Temp(i1) );
end
elseif Problem==3
  X1=load('nodes_2D_tri_linear');
X=X1(:,2:3);
T1=load('elem_2D_tri_linear');
T=T1(:,2:4);
nnode=length(X);
nelem=length(T);
for i=1:1:nnode
    X(i,3)=0;
end
	% printing heading to file
f=fopen('MyParaviewFile_2D_tri_linear.vtk','w');
fprintf(f,'# vtk DataFile Version 1.0\n');
fprintf(f,'ECM-CELL DIFFUSION-MECHANICS\n');
fprintf(f,'ASCII\n');
fprintf(f,'\n');
fprintf(f,'DATASET UNSTRUCTURED_GRID\n');
fprintf(f,'%s %8i %s\n','POINTS', nnode ,'float');

%%%%%%%%%%%%%%%%%%%%%% WRITING COORDINATES OF NODES %%%%%%%%%%%%%%%%%%%%%%%%%

	% printing coordinates
for i1=1:nnode
           fprintf(f,'%14.8E %14.8E %14.8E\n',X(i1,1:3));
end
fprintf(f,'\n');

%%%%%%%%%%%%%%%%%%%%%% WRITING CONNECTIVITY OF NODES %%%%%%%%%%%%%%%%%%%%%%%%%

fprintf(f,'%s %8i %8i\n','CELLS', nelem , nelem*4);
for i1=1:nelem
                
            fprintf(f,'%4i  %10i  %10i %10i\n',3,T(i1,1)-1,T(i1,2)-1,T(i1,3)-1);
end
fprintf(f,'\n');
fprintf(f,'%s %8i\n','CELL_TYPES', nelem);
for i1=1:nelem
            fprintf(f,' %4i ', 5);
end
fprintf(f,'\n');          

%%%%%%%%%%%%%%%%%%%%%% WRITING VARIABLES %%%%%%%%%%%%%%%%%%%%%%%%%
fprintf(f,'%s %8i\n','POINT_DATA', nnode);
fprintf(f,'SCALARS Ux float\n');
fprintf(f,'LOOKUP_TABLE default\n');
for i1=1:nnode
           fprintf(f,'%14.8E\n',Temp(i1) );
end

elseif Problem==4
    X1=load('nodes_2D_tri_quad');
X=X1(:,2:3);
T1=load('elem_2D_tri_quad');
T=T1(:,2:7);
nnode=length(X);
nelem=length(T);
for i=1:1:nnode
    X(i,3)=0;
end
	% printing heading to file
f=fopen('MyParaviewFile_2D_tri_quad.vtk','w');
fprintf(f,'# vtk DataFile Version 1.0\n');
fprintf(f,'ECM-CELL DIFFUSION-MECHANICS\n');
fprintf(f,'ASCII\n');
fprintf(f,'\n');
fprintf(f,'DATASET UNSTRUCTURED_GRID\n');
fprintf(f,'%s %8i %s\n','POINTS', nnode ,'float');

%%%%%%%%%%%%%%%%%%%%%% WRITING COORDINATES OF NODES %%%%%%%%%%%%%%%%%%%%%%%%%

	% printing coordinates
for i1=1:nnode
           fprintf(f,'%14.8E %14.8E %14.8E\n',X(i1,1:3));
end
fprintf(f,'\n');

%%%%%%%%%%%%%%%%%%%%%% WRITING CONNECTIVITY OF NODES %%%%%%%%%%%%%%%%%%%%%%%%%

fprintf(f,'%s %8i %8i\n','CELLS', nelem , nelem*7);
for i1=1:nelem
                
            fprintf(f,'%4i  %10i  %10i %10i %10i %10i %10i\n',6,T(i1,1)-1,T(i1,2)-1,T(i1,3)-1,T(i1,4)-1,T(i1,5)-1,T(i1,6)-1);
end 
fprintf(f,'\n');
fprintf(f,'%s %8i\n','CELL_TYPES', nelem);
for i1=1:nelem
            fprintf(f,' %4i ', 22);
end
fprintf(f,'\n');          

%%%%%%%%%%%%%%%%%%%%%% WRITING VARIABLES %%%%%%%%%%%%%%%%%%%%%%%%%
fprintf(f,'%s %8i\n','POINT_DATA', nnode);
fprintf(f,'SCALARS Ux float\n');
fprintf(f,'LOOKUP_TABLE default\n');
for i1=1:nnode
           fprintf(f,'%14.8E\n',Temp(i1) );
end

elseif Problem==5
  X1=load('nodes_3D_quad_lin.txt');
X=X1(:,2:4);
T1=load('elem_3D_quad_lin.txt');
T=T1(:,2:9);
nnode=length(X);
nelem=length(T);
% for i=1:1:nnode
%     X(i,3)=0;
% end
	% printing heading to file
f=fopen('MyParaviewFile_3D_quad_lin.vtk','w');
fprintf(f,'# vtk DataFile Version 1.0\n');
fprintf(f,'ECM-CELL DIFFUSION-MECHANICS\n');
fprintf(f,'ASCII\n');
fprintf(f,'\n');
fprintf(f,'DATASET UNSTRUCTURED_GRID\n');
fprintf(f,'%s %8i %s\n','POINTS', nnode ,'float');

%%%%%%%%%%%%%%%%%%%%%% WRITING COORDINATES OF NODES %%%%%%%%%%%%%%%%%%%%%%%%%

	% printing coordinates
for i1=1:nnode
       fprintf(f,'%14.8E %14.8E %14.8E\n',X(i1,1:end));
end
fprintf(f,'\n');

%%%%%%%%%%%%%%%%%%%%%% WRITING CONNECTIVITY OF NODES %%%%%%%%%%%%%%%%%%%%%%%%%

fprintf(f,'%s %8i %8i\n','CELLS', nelem , nelem*9);
for i1=1:nelem
                
            fprintf(f,'%4i %10i  %10i  %10i  %10i  %10i %10i %10i %10i\n',8,T(i1,1)-1,T(i1,2)-1,T(i1,3)-1,T(i1,4)-1,T(i1,5)-1,T(i1,6)-1,T(i1,7)-1,T(i1,8)-1);
end
fprintf(f,'\n');
fprintf(f,'%s %8i\n','CELL_TYPES', nelem);
for i1=1:nelem
            fprintf(f,' %4i ', 12);
end
fprintf(f,'\n');          

%%%%%%%%%%%%%%%%%%%%%% WRITING VARIABLES %%%%%%%%%%%%%%%%%%%%%%%%%
fprintf(f,'%s %8i\n','POINT_DATA', nnode);
fprintf(f,'SCALARS Ux float\n');
fprintf(f,'LOOKUP_TABLE default\n');
for i1=1:nnode
           fprintf(f,'%14.8E\n',Temp(i1) );
end

elseif Problem==6
    X1=load('nodes_3D_quad_quad.txt');
X=X1(:,2:4);
T1=load('elem_3D_quad_quad.txt');
T=T1(:,2:21);
nnode=length(X);
nelem=length(T);
% for i=1:1:nnode
%     X(i,3)=0;
% end
	% printing heading to file
f=fopen('MyParaviewFile_3D_quad_quad.vtk','w');
fprintf(f,'# vtk DataFile Version 1.0\n');
fprintf(f,'ECM-CELL DIFFUSION-MECHANICS\n');
fprintf(f,'ASCII\n');
fprintf(f,'\n');
fprintf(f,'DATASET UNSTRUCTURED_GRID\n');
fprintf(f,'%s %8i %s\n','POINTS', nnode ,'float');

%%%%%%%%%%%%%%%%%%%%%% WRITING COORDINATES OF NODES %%%%%%%%%%%%%%%%%%%%%%%%%

	% printing coordinates
for i1=1:nnode
       fprintf(f,'%14.8E %14.8E %14.8E\n',X(i1,1:end));
end
fprintf(f,'\n');

%%%%%%%%%%%%%%%%%%%%%% WRITING CONNECTIVITY OF NODES %%%%%%%%%%%%%%%%%%%%%%%%%

fprintf(f,'%s %8i %8i\n','CELLS', nelem , nelem*21);
for i1=1:nelem
                
            fprintf(f,'%4i %10i  %10i  %10i  %10i  %10i  %10i  %10i  %10i  %10i  %10i  %10i  %10i  %10i  %10i  %10i  %10i  %10i %10i %10i %10i\n',20,T(i1,1)-1,T(i1,2)-1,T(i1,3)-1,T(i1,4)-1,T(i1,5)-1,T(i1,6)-1,T(i1,7)-1,T(i1,8)-1,T(i1,9)-1,T(i1,10)-1,T(i1,11)-1,T(i1,12)-1,T(i1,13)-1,T(i1,14)-1,T(i1,15)-1,T(i1,16)-1,T(i1,17)-1,T(i1,18)-1,T(i1,19)-1,T(i1,20)-1);
end
fprintf(f,'\n');
fprintf(f,'%s %8i\n','CELL_TYPES', nelem);
for i1=1:nelem
            fprintf(f,' %4i ', 25);
end
fprintf(f,'\n');          

%%%%%%%%%%%%%%%%%%%%%% WRITING VARIABLES %%%%%%%%%%%%%%%%%%%%%%%%%
fprintf(f,'%s %8i\n','POINT_DATA', nnode);
fprintf(f,'SCALARS Ux float\n');
fprintf(f,'LOOKUP_TABLE default\n');
for i1=1:nnode
           fprintf(f,'%14.8E\n',Temp(i1) );
end

elseif Problem==7
X1=load('nodes_3D_tri_lin.inp');
X=X1(:,2:4);
T1=load('elem_3D_tri_lin.txt');
T=T1(:,2:5);
nnode=length(X);
nelem=length(T);

	% printing heading to file
f=fopen('MyParaviewFile_3D_tri_linear.vtk','w');
fprintf(f,'# vtk DataFile Version 1.0\n');
fprintf(f,'ECM-CELL DIFFUSION-MECHANICS\n');
fprintf(f,'ASCII\n');
fprintf(f,'\n');
fprintf(f,'DATASET UNSTRUCTURED_GRID\n');
fprintf(f,'%s %8i %s\n','POINTS', nnode ,'float');

%%%%%%%%%%%%%%%%%%%%%% WRITING COORDINATES OF NODES %%%%%%%%%%%%%%%%%%%%%%%%%

	% printing coordinates
for i1=1:nnode
       fprintf(f,'%14.8E %14.8E %14.8E\n',X(i1,1:end));
end
fprintf(f,'\n');

%%%%%%%%%%%%%%%%%%%%%% WRITING CONNECTIVITY OF NODES %%%%%%%%%%%%%%%%%%%%%%%%%

fprintf(f,'%s %8i %8i\n','CELLS', nelem , nelem*5);
for i1=1:nelem
                
            fprintf(f,'%4i  %10i  %10i %10i %10i\n',4,T(i1,1)-1,T(i1,2)-1,T(i1,3)-1,T(i1,4)-1);
end
fprintf(f,'\n');
fprintf(f,'%s %8i\n','CELL_TYPES', nelem);
for i1=1:nelem
            fprintf(f,' %4i ', 10);
end
fprintf(f,'\n');          

%%%%%%%%%%%%%%%%%%%%%% WRITING VARIABLES %%%%%%%%%%%%%%%%%%%%%%%%%
fprintf(f,'%s %8i\n','POINT_DATA', nnode);
fprintf(f,'SCALARS Ux float\n');
fprintf(f,'LOOKUP_TABLE default\n');
for i1=1:nnode
           fprintf(f,'%14.8E\n',Temp(i1) );
end

elseif Problem==8
X1=load('nodes_3D_tri_quad.txt');
X=X1(:,2:4);
T1=load('elem_3D_tri_quad.inp');
T=T1(:,2:11);
nnode=length(X);
nelem=length(T);
% for i=1:1:nnode
%     X(i,3)=0;
% end
	% printing heading to file
f=fopen('MyParaviewFile_3D_tri_quad.vtk','w');
fprintf(f,'# vtk DataFile Version 1.0\n');
fprintf(f,'ECM-CELL DIFFUSION-MECHANICS\n');
fprintf(f,'ASCII\n');
fprintf(f,'\n');
fprintf(f,'DATASET UNSTRUCTURED_GRID\n');
fprintf(f,'%s %8i %s\n','POINTS', nnode ,'float');

%%%%%%%%%%%%%%%%%%%%%% WRITING COORDINATES OF NODES %%%%%%%%%%%%%%%%%%%%%%%%%

	% printing coordinates
for i1=1:nnode
       fprintf(f,'%14.8E %14.8E %14.8E\n',X(i1,1:end));
end
fprintf(f,'\n');

%%%%%%%%%%%%%%%%%%%%%% WRITING CONNECTIVITY OF NODES %%%%%%%%%%%%%%%%%%%%%%%%%

fprintf(f,'%s %8i %8i\n','CELLS', nelem , nelem*11);
for i1=1:nelem
                
            fprintf(f,'%4i  %10i  %10i %10i  %10i %10i  %10i %10i  %10i %10i %10i\n',10,T(i1,1)-1,T(i1,2)-1,T(i1,3)-1,T(i1,4)-1,T(i1,5)-1,T(i1,6)-1,T(i1,7)-1,T(i1,8)-1,T(i1,9)-1,T(i1,10)-1);
end
fprintf(f,'\n');
fprintf(f,'%s %8i\n','CELL_TYPES', nelem);
for i1=1:nelem
            fprintf(f,' %4i ', 24);
end
fprintf(f,'\n');          

%%%%%%%%%%%%%%%%%%%%%% WRITING VARIABLES %%%%%%%%%%%%%%%%%%%%%%%%%
fprintf(f,'%s %8i\n','POINT_DATA', nnode);
fprintf(f,'SCALARS Ux float\n');
fprintf(f,'LOOKUP_TABLE default\n');
for i1=1:nnode
           fprintf(f,'%14.8E\n',Temp(i1) );
end

end

    toc
