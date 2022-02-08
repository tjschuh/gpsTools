function varargout=mat2mod(files)
% [dmat,tmax]=MAT2MOD(files)
%
% Given Precise Point Position time series of different units (e.g.,
% produced by PRD2MAT), makes them all start and end at the same time and
% inserts NaNs for times where no data were processed.
%
% INPUT:
% 
% files        cell with MAT-filename strings containing data structures
%
% OUTPUT:
%
% dmat         higher-dimensional structure with modified input structures
% tmax         two time strings with the inclusive range
%
% SEE ALSO:
%
% PRD2MAT
% 
% Originally written by tschuh-at-princeton.edu, 11/12/2021
% Last modified by tschuh-at-princeton.edu, 11/15/2021
% Last modified by fjsimons-at-alum.mit.edu, 02/08/2022

% Non-array variables to exclude from the tabling procedure
drem={'xyzunit','latlonunit','utmunit','heightunit','satlabels','utmzone'};
% Need to reconcile with earlier versions that mixed up the names
drem={'xyzunit','lonlatunit','utmunit','heightunit','satlabels','utmzone'};

for i=1:length(files)
    load(files{i});
    % Use (RE)TIME(TABLE) to fill in time skips/data gaps with NaNs
    [d,fnd]=retimes(d,drem,{'secondly','fillwithmissing'});
    % Assemble for later use
    dmat(i) = d;
end

% The longest common overlapping segment goes from the latest start to the
% earliest end of any file. Assuming you know there IS a field t...
for i=1:length(files)
   B(i)=dmat(i).t(1);
   E(i)=dmat(i).t(end);
end
% Earliest beginning, and latest end
tmax=[min(B) max(E)];

% Latest beginning, and earliest end
B=max(B);
E=min(E);

% Now select only the strictly interior overlapping points
for i=1:length(files)
  % Update the NON-TIME fields! Time is first, use it last
  for k=length(fnd):-1:1
    dmat(i).(fnd{k})=dmat(i).(fnd{k})([dmat(i).t>=B & dmat(i).t<=E],:);
  end
end

% Variable output
varns={dmat,tmax};
varargout=varns(1:nargout);
