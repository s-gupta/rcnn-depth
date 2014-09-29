function front=paretoGroup(X)

% PARETOGROUP  To get the Pareto Front from a given set of points.
% synopsis:           front =paretoGroup (objectiveMatrix)
% where:
%   objectiveMatrix: [number of points X number of objectives] array
%   front:           [number of points X 1] logical vector to indicate if ith
%                    point belongs to the Pareto Front (true) or not (false).
%
% by Yi Cao, Cranfield University, 31 June 2007
% 
% Identify the Pareto Front from a set of points in objective space is the 
% most important and also the most time-consuming task in multi-objective 
% optimization. This code splits the given objective set into several 
% smaller groups to be examined by the efficient paretofront algorithm. 
% Then, the Pareto Fronts of each group are combined as one set to be 
% checked by the paretofront algorithm to determine the overall Pareto
% Front. In this way, the overal computation time can be reduced about
% half.
%
% Example:
% X = rand(1000000,4);
% t0 = cputime;
% Y1=paretoGroup(X); %mex implementation without sorting.
% t1=cputime - t0;
% t0 = cputime;
% Y2=paretofront(X);
% t2=cputime - t0;
% isequal(Y1,Y2)    %shoudl be 1
% disp([t1 t2])     
% Computation time based on Intel(R) Core(TM)2 CPU T2500 @ 2.0GHz, 2.0 GB of RAM
% 0.6844    1.4404
%

[m,n]=size(X);
groupcut=floor(2^13/n);
gRoup=max(1,ceil(m/groupcut));
front=false(m,1);
for k=1:gRoup
    z0=(k-1)*groupcut;
    z=(z0+1):min(z0+groupcut,m);
    front(z)=paretofront(X(z,:));
end
if gRoup>1
    front(front)=paretofront(X(front,:));
end