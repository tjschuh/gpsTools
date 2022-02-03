function varargout=mat2mod(files)
% [dmat,tmax]=MAT2MOD(files)
%
% Given Precise Point Position time series of (four) different units, makes
% them all start and end at the same time and inserts NaNs for times where
% no data were processed.
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
% EXAMPLE:
%
% dmat=mat2mod(files);
% 
% Originally written by tschuh-at-princeton.edu, 11/12/2021
% Last modified by tschuh-at-princeton.edu, 11/15/2021
% Last modified by fjsimons-at-alum.mit.edu, 02/01/2022

% Non-array variables to exclude from the tabling procedure
drem={'xyzunit','lonlatunit','utmunit','heightunit','satlabels'};
% need to fix prd2mat.m for rtk: latlonunit --> lonlatunit
%drem={'xyzunit','latlonunit','utmunit','heightunit','satlabels'};

for i=1:length(files)
    load(files{i});
    % Only get the wanted variables, turn them into a timetable
    % Use (RE)TIME(TABLE) to fill in time skips/data gaps with NaNs
    tt=retime(table2timetable(struct2table(rmfield(d,drem))),'secondly','fillwithmissing');
    % Reassign to the old STRUCT
    d.t=tt.t;
    % These are the variables we kept
    fnd=fieldnames(rmfield(d,drem));
    % And they appear shifted by one in the time table
    fnt=fieldnames(tt);
    % So now you put them back in the right place as a struct
    for j=2:length(fnd)
      d.(fnd{j})=tt.(fnt{j-1});
    end    
    % Assemble for later use
    dmat(i) = d;
end

% About to change this to: all files start and end with a datum, but
% in-between their values, there could be skips, which now have been
% fixed by setting them to NaNs. The longest common overlapping segment
% goes from the latest start to the earliest end of any file.
% Assuming you know there IS a field t.
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
