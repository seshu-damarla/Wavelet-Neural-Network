%dvgpro%

% Function: model=HSPXY(X,Y,Ncal)

% Input:
% X, matrix of predictor variables.
% Y, vector of the response variable.
% Ncal, number of objects to be selected to model set.

% Output:
% model, vector of objects selected to model set.

function [model,test]=HSPXY(X,Y,Ncal)

dminmax = zeros(1,Ncal); % Initializes the vector of minimum distances.

M= size(X,1); % Number of rows in X.
samples = 1:M;

Dx = zeros(M, M);  % Initializes the matrix of X Euclidean distances.
Dy = zeros(M,M);   % Initializes the matrix of Y Euclidean distances.
Dc = zeros(M, M);  % Initializes the matrix of X adverse cosine distances.

for i = 1:M-1
    xa =X(i,:);
    ya = Y(i,:);
    for j = i + 1:M
        xb =X(j,:);
        yb = Y(j,:);     
        xab=[xa;xb];
        Dx(i, j) = norm(xa-xb);
        Dy(i, j) = norm(ya-yb);
        Dc(i,j)=abs(1-pdist(xab,'cosine'));
    end
end
Dxmax =max(max(Dx));
Dymax =max(max(Dy));
Dcmax=max(max(Dc));
%Combine the subtraction of cosine distance 
D=1+ Dx/Dxmax -Dc/Dcmax +Dy/Dymax;

% D is an upper triangular matrix.
% D(i,j) is the distance between objects i and j (j > i).

[maxD,index_row] =max(D);
% maxD is a row vector containing the largest element for each column of D.
% index row is the row in which the largest element of the column if found.
[dummy, index_column] =max(maxD);
% index_column is the column containing the largest element of matrix D.
model(1) = index_row(index_column);
model(2) = index_column;
for i = 3:Ncal
    pool = setdiff(samples,model);
    % Pool is the index set of the samples that have not been selected yet.
    dmin = zeros(1,M-i + 1);
    % dmin will store the minimum distance of each sample in 'pool' with respect to the previously selected samples.
    for j = 1:(M-i + 1)
        indexa = pool(j);
        d = zeros(1,i-1);
        for k = 1:(i-1)
            indexb =model(k);
            if indexa < indexb
                d(k) =D(indexa,indexb);
            else                
                d(k) =D(indexb,indexa);
            end            
        end       
        dmin(j) = min(d);
 
    end
    % At each iteration, the sample with the largest dmin value is selected.
   [dummy,index] =max(dmin);
   model(i) = pool(index);
end
if nargout==2
    test=samples;
    test(model)=[];
end

