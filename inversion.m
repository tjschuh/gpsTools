function [sol,tru] = inversion(unit1file,unit2file,unit3file,unit4file)
%
% use trilateration method (not my code) to calculate beacon
% location on seafloor from slant ranges and ship locations
% can also use RecTrilateration.m to see if it improves result
%
% INPUT:
%
% unit1file     mat file containing data collected by unit 1
% unit2file     mat file containing data collected by unit 2
% unit3file     mat file containing data collected by unit 3
% unit4file     mat file containing data collected by unit 4
%
% OUTPUT:
%
% sol           calculated [x y z] beacon location [m]
% tru           actual [x y z] beacon location [m]
%
% Originally written by tschuh-at-princeton.edu, 11/30/2021
% Last modified by tschuh-at-princeton.edu, 01/04/2022

% use mat2mod to convert data to all be same time spans with no time gaps
[d1,d2,d3,d4] = mat2mod(unit1file,unit2file,unit3file,unit4file);

% true beacon location
% eventually this wont be known
dogx = 1.979e6;
dogy = -5.074e6;
dogz = 3.30385e6; %~5 km below surface
tru = [dogx dogy dogz];

% find a location (x,y,z) for every second 
for i=1:length(d1.xyz)
    % GPS locations [m]
    % ignoring d4 for now b/c it's bad!
    P1(:,i) = d1.xyz(i,:)';
    P2(:,i) = d2.xyz(i,:)';
    P3(:,i) = d3.xyz(i,:)';

    % compute slant ranges [m] for d1,d2,d3
    % eventually need to use slant times and velocity
    % profile to find slant ranges (sr = v*st)
    sr1 = sqrt((d1.xyz(i,1) - dogx).^2 + (d1.xyz(i,2) - dogy).^2 + (d1.xyz(i,3) - dogz).^2);
    sr2 = sqrt((d2.xyz(i,1) - dogx).^2 + (d2.xyz(i,2) - dogy).^2 + (d2.xyz(i,3) - dogz).^2);
    sr3 = sqrt((d3.xyz(i,1) - dogx).^2 + (d3.xyz(i,2) - dogy).^2 + (d3.xyz(i,3) - dogz).^2);

    P = [P1(:,i) P2(:,i) P3(:,i)];
    S = [sr1 sr2 sr3];

    % if there are any NaN values, dont do Trilateration, it wont run
    if isnan(sr1)==1 | isnan(sr2)==1 | isnan(sr3)==1
        solmat(i,:) = [NaN NaN NaN];
    else
        % use Trilateration code from "An Algebraic Solution to the Multilateration Problem"
        % soltype = 1 or 2 (2 is more accurate)
        soltype = 2;
        sol = Trilateration(P,S,diag(ones(1,3)),soltype);
        solmat(i,:) = [sol(2,1) sol(3,1) sol(4,1)];
    end
end

% for final location take average of all locations (maybe not best way)
% could make histogram here
sol = real(nanmean(solmat,1));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function sol = Trilateration(P,S,W,soltype)

% np = # of stations
% ns = # of slant ranges
[mp,np] = size(P);
ns = length(S);
% np and ns must match
if (ns ~= np)
    error('Number of reference points and distances are different');
end

A=[]; b=[];
for i1=1:np
    x = P(1,i1); y = P(2,i1); z = P(3,i1);
    s = S(i1);
    A = [A ; 1 -2*x  -2*y  -2*z]; 
    b= [b ; s^2-x^2-y^2-z^2 ];
end

if np == 3
    % Gaussian elimination
    Xp= A\b;
    % or Xp=pinv(A)*b; 
    % the matrix  inv(A'*A)*A' or inv(A'*C*A)*A'*C or pinv(A)
    % depend only on the reference points
    % it could be computed only once
    xp = Xp(2:4,:);
    Z = null(A,'r');
    z = Z(2:4,:);
    if rank(A) == 3
        % Polynomial coefficients
        a2 = z(1)^2 + z(2)^2 + z(3)^2 ;
        a1 = 2*(z(1)*xp(1) + z(2)*xp(2) + z(3)*xp(3))-Z(1);
        a0 = xp(1)^2 +  xp(2)^2+  xp(3)^2-Xp(1);
        p = [a2 a1 a0];
        t = roots(p);

        % Solutions
        if soltype == 1
            sol = Xp + t(1)*Z;
        elseif soltype == 2    
            sol = Xp + t(2)*Z;
        end
    end

elseif np > 3
% Particular solution

    if W ~= diag(ones(1,length(W)))
        C = W'*W;
        Xpdw = inv(A'*C*A)*A'*C*b; % Solution with Weights Matrix
    else
        Xpdw = pinv(A)*b; % Solution without Weights Matrix
    end
 
    % the matrix  inv(A'*A)*A' or inv(A'*C*A)*A'*C or pinv(A)
    % depend only on the reference points
    % it could be computed only once
    if soltype == 1
        sol = Xpdw;
    elseif soltype == 2
        sol = Xpdw;
    end

else % np <= 3
    error('Not enough stations');
end