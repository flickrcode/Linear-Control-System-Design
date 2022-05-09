%% SSY285 Linear Control System Design:  Assignment - 1
%% Group 11 - Manikanta Venkatesh, Chintalapudi Adhitya Reddy and Fikri Farhan Witjaksono

clc
clear all
close all

%% Question d - Calculation of Eigen values of A
R = 1;
Ke = 0.1;
Kt = 0.1;
J1 = 1e-5;
J2 = 4e-5;
Bf = 2e-3;
D1 = 20;
D2 = 2;

A = [ 0 0 0 1 0 ;
     0 0 0 0  1;
     0 (D2/Bf) (-D2/Bf) 0 0 ;
     (-D1/J1) (D1/J1) 0 (-(Kt*Ke)/(R*J1)) 0;
     (D1/J2) ((-D1-D2)/J2) (D2/J2) 0 0];
 
Eigenvalue = eig(A)

%  Calculation of poles and zeroes for 1st case
B1= [0 0;
    0 0;
    0 (1/Bf);
   (Kt/(R*J1)) 0 ;
     0 0];

C1 = [0 1 0 0 0;
     0 0 0 0 1];
  
D1 = zeros(2);

G1 = ss(A,B1,C1,D1);
s1 = tf(G1);

for n=1:1:2
    for i=1:1:2
[zz1,pp1,kk1] = tf2zp(s1.numerator{n,i},s1.denominator{n,i})
    end
end

% Calculation of poles and zeroes for 2nd case
B2= [0 0;
    0 0;
    0 (1/Bf);
   (Kt/(R*J1)) 0 ;
     0 0];

C2 = [ 0 0 0 ((-Ke)/R) 0;
     0 (D2/Bf) (-D2/Bf) 0 0];
  
D2 = [(1/R) 0;
    0 (1/Bf)];

G2 = ss(A,B2,C2,D2);
s2 = tf(G2);

for n=1:1:2
    for i=1:1:2
[zz2,pp2,kk2] = tf2zp(s2.numerator{n,i},s2.denominator{n,i})
    end
end

%% Question f - Calculation of poles and transmission zeroes with no external torque
R = 1;
Ke = 0.1;
Kt = 0.1;
J1 = 1e-5;
J2 = 4e-5;
Bf = 2e-3;
D1 = 20;
D2 = 2;

A = [ 0 0 0 1 0 ;
     0 0 0 0  1;
     0 (D2/Bf) (-D2/Bf) 0 0 ;
     (-D1/J1) (D1/J1) 0 (-(Kt*Ke)/(R*J1)) 0;
     (D1/J2) ((-D1-D2)/J2) (D2/J2) 0 0];
 
B3= [0;
    0;
    0 ;
   (Kt/(R*J1)) ;
     0];

C3 = [ 0 0 0 ((-Ke)/R) 0;
     0 (D2/Bf) (-D2/Bf) 0 0];
  
D3 = [(1/R); 0];

G3 = ss(A,B3,C3,D3);
 s3 = tf(G3);

for n=1:1:2
    [zz3,pp3,kk3] = tf2zp(s3.numerator{n,1},s3.denominator{n,1})
end

