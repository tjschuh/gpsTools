function varargout=retimei(numval,jmp)
% [nnumval,numval]=RETIMEI(numval,jmp)
%
% An index version of RETIME for a two-column data set
%
% INPUT:
%
% numval    A Mx2 matrix with index-value pairs
% jmp       If the index (column one) exceeds this value you insert NaN
%
% OUTPUT:
%
% nnumval    A new array with the missing indices replaced
% numval     The thing you had put in
%
% EXAMPLE:
%
% retimei([],2)
%
% Last modified by fjsimons-at-alum.mit.edu, 02/08/2022

defval('numval',[[1 2 3 6 7 8 13 15 17 18 19]',[1 2 3 4 5 6 7 8 9 10 11]'])
defval('jmp',2);

% Just play with it...
dtg=diff(numval(:,1));
% Find where the big jumps in the independent variable are
bigj=find(dtg>jmp);
% Insert the right amount of NaNs to take care of those jumps
news1=insert(numval(:,1),NaN,gamini(bigj+1,dtg(bigj)-1))';

% news2=insert(numval(:,2),NaN,gamini(bigj+1,dtg(bigj)-1))';
% Then either redo the line above or be smarter
news2=nan(size(news1));
news2(~isnan(news1))=numval(:,2);

% The newfangled thing
nnumval=[news1 news2];

% Optional output
varns={nnumval,numval};
varargout=varns(1:nargout);
