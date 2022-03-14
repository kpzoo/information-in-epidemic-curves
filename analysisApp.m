% Examine various determinants to show extension of results in the Appendix
% for section "Information if case source times are never reported"
clearvars; clc; close all; 

% Assumptions and notes
% - assumes delays cannot be separated and builds up to 4x4 FIM
% - X1-X4 are the algebraic matices of relevance to the Appendix
% - here A = R1 so a1/A terms is same as (1,1) in the manuscript

% Coefficients of R terms
syms a1 b1 b2 c1 c2 c3 d1 d2 d3 d4
% R terms
syms A B C D

% Except for first case all should be identifiable
I2 = [b1^2 b1*b2; b1*b2 b2^2];
I3 = [c1^2 c1*c2 c1*c3; c1*c2 c2^2 c2*c3; c1*c3 c2*c3 c3^2];
I4 = [d1^2 d1*d2 d1*d3 d1*d4; d1*d2 d2^2 d2*d3 d2*d4; d1*d3 d2*d3 d3^2 d3*d4; d1*d4 d2*d4 d3*d4 d4^2];

% Get determinants proving non-identifiability
D2 = eval(det(I2)); D3 = eval(det(I3)); D4 = eval(det(I4));
disp(['Determinants are: ' num2str([D2 D3 D4])]);

% Sum up in stages - combines information along epidemic
J1 = a1/A;
J2 = I2/B; J2(1,1) = J2(1,1) + J1;
J3 = I3/C; J3(1:2, 1:2) = J3(1:2, 1:2) + J2;
J4 = I4/D; J4(1:3, 1:3) = J4(1:3, 1:3) + J3;

% Display all determinants
Ds = [det(J1) det(J2) det(J3) det(J4)];
pretty(Ds)

% Some simple manipulation to expose determinant
X2 = J2; 
X2(:,1) = X2(:,1) - (b1/b2)*X2(:,2);
X3 = J3; 
X3(:,1) = X3(:,1) - (c1/c3)*X3(:,3);
X3(:,2) = X3(:,2) - (c2/c3)*X3(:,3); 
X3(:,1) = X3(:,1) - (b1/b2)*X3(:,2);
X4 = J4;
% Remove all d parts
X4(:,1) = X4(:,1) - (d1/d4)*X4(:,4);
X4(:,2) = X4(:,2) - (d2/d4)*X4(:,4);
X4(:,3) = X4(:,3) - (d3/d4)*X4(:,4);
% Remove all c parts as in X3
X4(:,1) = X4(:,1) - (c1/c3)*X4(:,3);
X4(:,2) = X4(:,2) - (c2/c3)*X4(:,3); 
% Remove all b parts as in X2
X4(:,1) = X4(:,1) - (b1/b2)*X4(:,2);

% Show precise comparison with Appendix for 3x3
syms alpha01 alpha02 alpha03 alpha11 alpha12 alpha21 beta1 beta2 beta3
A = beta1/alpha01; B = beta2; C = beta3;
a1 = alpha01; b2 = alpha02; c3 = alpha03;
b1 = alpha11; c2 = alpha12; c1 = alpha21;
X3app = simplify(eval(X3)); pretty(X3app)