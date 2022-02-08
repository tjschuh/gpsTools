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

% Save to MAT

% Make a plot if there is no output? Watch out for implicit vector sums
%plot([tags(:,1)-tags(1,1)]*1e9+tags(:,2)-[1:size(tags(:,1),1)]')
%plot(tags(:,1)*1e9+tags(:,2))

% Optional output
varns={tags};
varargout=varns(1:nargout);

d e f i n i t e l y n o t r e a d y y e t

% 1 and 3 have the most points
load DOG1-camp
[sr]=gps2rng({'Unit3-camp.mat'});
inta=floor(sr/1500);
plot(sr/1600)
hold on
ofs=0
for ofs=000:3000
  p=plot(tags(:,2)/1e9+inta([1:size(tags,1)]+ofs),'k');
  title(num2str(ofs))                              
  pause
  delete(p)
end
