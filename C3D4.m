function [n,w,xi,N,dNdxi]=C2D4
%====================== No. integration points =============================
%
%   Defines the number of integration points:be used for
%   each element type
%

n = 1;
ncoord=3;  
nodes=4;
npoints = 1;
%
%====================== INTEGRATION POINTS ==================================
%
%   Defines positions of integration points
%
xi = zeros(ncoord,npoints);
         xi(1,1) = 0.25;
         xi(2,1) = 0.25;
         xi(3,1) = 0.25;

         w = [1./6];
%
%================= SHAPE FUNCTIONS ==================================
%
%        Nij: Shape functions of the Int Point i [4x4] Ni [4x1]

N = zeros(n,nodes);
xi = xi';
for i1=1:n
       N(i1,1) = xi(i1,1);
       N(i1,2) = xi(i1,2);
       N(i1,3) = xi(i1,3);
       N(i1,4) = 1.-xi(i1,1)-xi(i1,2)-xi(i1,3);
end
%
%================= SHAPE FUNCTION DERIVATIVES ======================
%
%        Nij,r: Dev of shape functions of the Int Point i [4x8]
%        [2*i-1 2*i] => dNi,r [4x2]
dNdxi = zeros(ncoord*n,nodes);
for i1=1:n
       dNdxi(1,1) = 1.;
       dNdxi(2,2) = 1.;
       dNdxi(3,3) = 1.;
       dNdxi(4,1) = -1.;
       dNdxi(4,2) = -1.;
       dNdxi(4,3) = -1.;
end
end
%
