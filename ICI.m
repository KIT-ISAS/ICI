function [c,C,omega] = ICI(xA,CA,xB,CB)
% Inverse Covariance Intersection 
%
% This function implements the ICI algorithm 
% and fuses two estimates (xA,CA) and (xB,CB). 
% It provides the fusion result (c,C) and the 
% value of omega, which minimizes the trace of C. 

f = @(w)trace(inv(inv(CA)+...
        inv(CB)-inv(w*CA+(1-w)*CB)));

omega = fminbnd(f,0,1,optimset('Display','off'));
C = inv(inv(CA)+inv(CB)-...
        inv(omega*CA+(1-omega)*CB));

% computations of the gains
KICI = C*(inv(CA)-...
        omega*inv(omega*CA+(1-omega)*CB));
LICI = C*(inv(CB)-...
        (1-omega)*inv(omega*CA+(1-omega)*CB));

c = KICI*xA + LICI*xB;

