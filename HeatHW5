%% 
% Cole Stoumbaugh
clear;clc
% Given Times
t = [2.5 25 250] * 60;
% Given Lengths
L1 = 0.6;
L2 = 0.3;
Lc = 0.15;
% Other Given Properties
k = 30;
p = 7900;
Cp = 640;
h1 = 10;
h2 = 100;
% Properties from table 5.1
z1 = [0.3111,3.1731,6.2991,9.4354];
z2 = [0.8603,3.4256,6.4373,9.5293];

% Calculating Bi and Fo and T
Fo = (k * t)/(p*Cp*(L2^2));
Bi1 = (h1*L2)/k;
Bi2 = (h2*L2)/k;
T1 = (p*Cp*L2)/h1;
T2 = (p*Cp*L2)/h2;

%% Calculating theta-star
% Exact Solution
n = 1;
for m = 1:3
for n = 1:4
C(1,n) = (4*sin(z1(1,n)))/((2*z1(1,n)) + sin(2*z1(1,n)));
C(1,n+4) = (4*sin(z2(1,n)))/((2*z2(1,n)) + sin(2*z2(1,n)));
O11(m,n) = C(1,n)*exp(-(z1(1,n)^2)*Fo(1,m))*cos(z1(1,n));
O12(m,n) = C(1,n+4)*exp(-(z2(1,n)^2)*Fo(1,m))*cos(z2(1,n));
n = n+1;
end
O111(1,m) = sum(O11(m,:));
O112(1,m) = sum(O12(m,:));
m = m+1;
end

% L.C.
O21 = exp(-t/T1);
O22 = exp(-t/T2);

% First Term
O311 = C(1,1)*exp(-(z1(1,1)^2)*Fo)*cos(z1(1,1));
O312 = C(1,5)*exp(-(z2(1,1)^2)*Fo)*cos(z2(1,1));

% Semi-infinite Solid
for y = 1:3
O41(1,y) = exp((Bi1^2)*Fo(1,y))*erfc(Bi1*sqrt(Fo(1,y)));
O42(1,y) = exp((Bi2^2)*Fo(1,y))*erfc(Bi2*sqrt(Fo(1,y)));
y = y+1;
end
O = [O111, O112; O21, O22; O311, O312; O41, O42];

%% Table

FoNum = ["Exact Soultion"; "Lump Capacitance"; "First Term"; "Semi-infinite"];
Bi01Fo1 = O(:,1);
Bi01Fo2 = O(:,2);
Bi01Fo3 = O(:,3);
Bi02Fo1 = O(:,4);
Bi02Fo2 = O(:,5);
Bi02Fo3 = O(:,6);
ThetaStar = table(FoNum,Bi01Fo1,Bi01Fo2,Bi01Fo3,Bi02Fo1,Bi02Fo2,Bi02Fo3)
