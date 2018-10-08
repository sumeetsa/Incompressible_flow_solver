clc
clear all
X = importdata('xgrid.dat');
Y = importdata('ygrid.dat');
P = importdata('pressure.dat');
U = importdata('vel_u.dat');
V = importdata('vel_v.dat');
figure(1)
contourf(P)
figure(2)
contourf(U)
figure(3)
contourf(V)
figure(4)
streamslice(X,Y,U,V)