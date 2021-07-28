function b = myMontage(a,rc)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% myMontage.m - Richard Dortch
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Taken from MATLAB's montage.m.  Allows for user specified
% or calculated # of rows and colums for montage and returns 
% it.
%
% Input: a - 3D array with images to be used in montage
%        rc(optional) - [row column] of calculated montage
% Output: b - calculated montage
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Size of a
[nRows, nCols, nFrames] = size(a);

% Set number of row and colums in montage
if exist('rc') == 0 
    
% Columns to Rows to be one (square montage)
aspectRatio = 1;
nMontageCols = sqrt(aspectRatio * nRows * nFrames / nCols);

% Make sure montage rows and columns are integers
nMontageCols = ceil(nMontageCols);
nMontageRows = ceil(nFrames / nMontageCols);
else
    nMontageRows = rc(1);
    nMontageCols = rc(2);
end

% Create the montage image.
b = zeros(nMontageRows*nRows, nMontageCols*nCols);
rows = 1 : nRows;
cols = 1 : nCols;
for i = 0:nMontageRows-1
    for j = 0:nMontageCols-1,
        k = j + i * nMontageCols + 1;
        if k <= nFrames
            b(rows + i * nRows, cols + j * nCols, :) = a(:,:,k);
        else
            break
        end
    end
end
