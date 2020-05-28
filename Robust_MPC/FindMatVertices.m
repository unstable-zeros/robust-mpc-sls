function [matVertices] = FindMatVertices(dim, eps)
%FINDMATVERTICES find the vertices of an inf-norm ball of radius eps in the
%matrix space whose dimension is specified by dim.

if nargin < 2
    eps = 1;
end

% n: row size; m: column size.
n = dim(1); m = dim(2);

if eps == 0
   matVertices = cell(1,1);
   matVertices{1} = zeros(n, m);
   return;
end

if n == 1
   matVertices = cell(1, 2*m);
   for j = 1:m
       temp = zeros(1, m); temp(j) = eps;
       matVertices{2*j - 1} = temp;
       temp = zeros(1, m); temp(j) = -eps;
       matVertices{2*j} = temp;
   end
else
    lowdim = n - 1;
    newDim = [lowdim m];
    matVerticesLowDim = FindMatVertices(newDim, eps);
    
    numLowDimVertices = size(matVerticesLowDim, 2);
    curMatVertSize = 2*m*numLowDimVertices;
    matVertices = cell(1, curMatVertSize);
    for ii = 1:m
        for jj = 1:numLowDimVertices
            temp = zeros(1, m); temp(ii) = eps;
            matVertices{(ii-1)*2*numLowDimVertices + 2*jj - 1} = [temp; matVerticesLowDim{jj}];
            temp = zeros(1, m); temp(ii) = -eps;
            matVertices{(ii-1)*2*numLowDimVertices + 2*jj} = [temp; matVerticesLowDim{jj}];
        end
    end
        
end

end


