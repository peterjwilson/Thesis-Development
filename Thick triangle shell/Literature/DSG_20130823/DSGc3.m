clear all;
clc;

syms x1 x2 x3 y1 y2 y3 x y W1 W2 W3 PX1 PX2 PX3 PY1 PY2 PY3

% x1 = 0
% y1 = 0
% x2 = 1
% y2 = 0
% x3 = 0
% y3 = 1

A = [-1 -x1 -y1 -(x1*x1) -0.5*x1*y1 -0.5*x1*y1 -(y1*y1) x1*(1-y1) y1*(1-x1);
    -1 -x2 -y2 -(x2*x2) -0.5*x2*y2 -0.5*x2*y2 -(y2*y2) x2*(1-y2) y2*(1-x2);
    -1 -x3 -y3 -(x3*x3) -0.5*x3*y3 -0.5*x3*y3 -(y3*y3) x3*(1-y3) y3*(1-x3);
    0 1 0 2*x1 y1 0 0 0 0;
    0 1 0 2*x2 y2 0 0 0 0;
    0 1 0 2*x3 y3 0 0 0 0;
    0 0 1 0 0 x1 2*y1 0 0;
    0 0 1 0 0 x2 2*y2 0 0;
    0 0 1 0 0 x3 2*y3 0 0]

%Ainv = inv(A);

U = [W1; W2; W3; PX1; PX2; PX3; PY1; PY2; PY3]

AinvU = A\U;

simplify(AinvU)