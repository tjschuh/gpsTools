function varargout=dog2mat(diro,fname,xver)
% tags=DOG2MAT(diro,fname,xver)
%
% Reads in all the *.txt files with data from the geodesy board
% Note that standard Unix/Linux "more", "cat", "tail" etc do not work on
% these files but "vi" and "emacs" and MATLAB load handle them properly.
%
% INPUT
%
% diro        Some directory
% fname       The save file name
% xver        1 Makes various plots
% 
% OUTPUT:
%
% tags        The tags in two-column format
%
% SEE ALSO:
%
% DOG2RNG which performs the unwrapping.
%
% Last modified by fjsimons-at-alum.mit.edu, 08/14/2023

% Where data are being kep
defval('diro','/data1/seafloorgeodesy/GeodesyBoard/DOG1/camp/txt')

% Make the target filename
defval('fname','DOG1-camp.mat')

if exist(fname)==2
  load(fname)
else
  % Finds all the files, e.g. File037.txt, etc
  files=ls2cell(diro);
  
  % Read in all the files, which are comma-separated ASCII files
  tags=[];
  for index=1:length(files)
    d=load(files{index});
    tags=[tags ; d];
  end

  if xver==1
    % Make a plot if there is no output? Watch out for implicit vector sums
  end

  % Save the filename
  save(fname,'tags')
end

% Optional output
varns={tags};
varargout=varns(1:nargout);
  
