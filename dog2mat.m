function tags=dog2mat(diro)
% tags=DOG2MAT(diro)
%
% Reads in all the *.txt files with data from the geodesy board
%
% INPUT
%
% diro        Some directory
% 
% OUTPUT:
%
% tags        The tags in two-column format
%
% EXAMPLE:
%
% Ran this, then saved the outputs in a *mat
%
% Last modified by fjsimons-at-alum.mit.edu, 02/07/2022

% Finds all the files
files=ls2cell(diro);

tags=[];
for index=1:length(files)
  d=load(files{index});
  tags=[tags ; d];
end

% Optional output
varns={tags};
varargout=varns(1:nargout);

