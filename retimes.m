function d=retimes(d,drem,options)
% d=RETIMES(d,drem,options)
%
% Allows RETIME to work on a structure provided you specify what non-numeric
% fields to ignore. Output is again a struct.
%
% INPUT:
%
% d         A struct
% drem      A cell name with the fieldname strings that you want to ignore
% options   The arguments you want to pass onto RETIME (two for now)
%
% OUTPUT:
%
% d      The new struct
%
% Last modified by fjsimons-at-alum.mit.edu, 02/06/2022

defval('drem',{'xyzunit','lonlatunit','utmunit','heightunit','satlabels'});
defval('options',{'secondly','fillwithmissing'})

% Only get the wanted variables, turn them into a timetable
tt=retime(table2timetable(struct2table(rmfield(d,drem))),options{1},options{2});
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
