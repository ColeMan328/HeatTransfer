% Cole Stoumbaugh
clear;clc

%% Variable Setup
Alum = readmatrix("Lab3Data_aluminum_0.5.txt");
Ny1 = readmatrix("Lab3Data_nylon_0.5.txt");
Ny2 = readmatrix("Lab3Data_nylon_0.75.txt");
Ny3 = readmatrix("Lab3Data_nylon_1.0.txt");
Steel = readmatrix("Lab3Data_steel_0.5.txt");

tAlum = Alum(:,1);
tinfAlum = Alum(:,2)+273;
dtinfAlum = Alum(:,3);
t1Alum = Alum(:,4)+273;
dt1Alum = Alum(:,5);
tstAlum = Alum(1,4)+273;
D1Alum = 0.5/39.37;
pAlum = 2790;
cAlum = 880;
kAlum = 134;
Lc = (1/6)*D1Alum;

tNy1 = Ny1(:,1);
tinfNy1 = Ny1(:,2)+273;
dtinfNy1 = Ny1(:,3);
t1Ny1 = Ny1(:,4)+273;
dt1Ny1 = Ny1(:,5);
tstNy1 = Ny1(1,4)+273;
D1Ny1 = 0.5/39.37;
pNy1 = 1140;
cNy1 = 1500;
kNy1 = 0.2;
LcNy1 = (1/6)*D1Ny1;

tNy2 = Ny2(:,1);
tinfNy2 = Ny2(:,2)+273;
dtinfNy2 = Ny2(:,3);
t1Ny2 = Ny2(:,4)+273;
dt1Ny2 = Ny2(:,5);
tstNy2 = Ny2(1,4)+273;
D1Ny2 = 0.75/39.37;
pNy2 = 1140;
cNy2 = 1500;
kNy2 = 0.2;
LcNy2 = (1/6)*D1Ny2;

tNy3 = Ny3(:,1);
tinfNy3 = Ny3(:,2)+273;
dtinfNy3 = Ny3(:,3);
t1Ny3 = Ny3(:,4)+273;
dt1Ny3 = Ny3(:,5);
tstNy3 = Ny3(1,4)+273;
D1Ny3 = 1/39.37;
pNy3 = 1140;
cNy3 = 1500;
kNy3 = 0.2;
LcNy3 = (1/6)*D1Ny3;

tSteel = Steel(:,1);
tinfSteel = Steel(:,2)+273;
dtinfSteel = Steel(:,3);
t1Steel = Steel(:,4)+273;
dt1Steel = Steel(:,5);
tstSteel = Steel(1,4)+273;
D1Alum = 0.5/39.37;
pSteel = 7870;
cSteel = 486;
kSteel = 51.9;

x = 0.367;
L = D1Alum;
o = 5.67*10^(-8);
e = 0.8;
B = 3.12*10^(-3);
v = 1.807*10^(-5);
u = 1.962*10^(-5);
c = 1.0063 * 10^3;
k = 2.816 * 10 ^-2;

%% Experimental Alum
n = 1;
m = 1;
while n <= numel(t1Alum)

O1Alum(n,1) = (t1Alum(n,1) - tinfAlum(n,1));
O2Alum(n,1) = (tstAlum - tinfAlum(n,1));
OAlum(n,1) = (O1Alum(n,1))/(O2Alum(n,1));
if (x+.001) > OAlum(n,1) > (x-.001)
    TouAlum(m,1) = tAlum(n,1);
    m = m+1;
end
n = n+1;
end
TouAlum = TouAlum(1,1);
hAlum = (pAlum*cAlum*Lc)/TouAlum;
BiAlum = (Lc*hAlum)/kAlum;

%% Theoretical Alum
tinfAlum = sum(tinfAlum)/numel(tinfAlum);
GrAlum = (9.81*B*(L^3)*(tstAlum-tinfAlum))/(v^2);
PrAlum = (c*u)/k;
RaAlum = PrAlum * GrAlum;
hconvAlum = (k/L)*(2+(0.589*RaAlum^(1/4))/(1+(0.469/PrAlum)^(9/16))^(4/9));
hradAlum = o*0.6*(tstAlum+tinfAlum)*((tstAlum^2)+(tinfAlum^2));
htotAlum = hconvAlum+hradAlum;
TouhAlum = (pAlum*cAlum*Lc)/htotAlum;
BihAlum = (htotAlum*Lc)/kAlum;
OhAlum = exp(-tAlum/TouhAlum);

perdiffTouAlum = ((TouAlum - TouhAlum)/TouhAlum)*100; 
perdiffBiAlum = ((BiAlum - BihAlum)/BihAlum)*100;

%% Experimental Ny1
n = 1;
m = 1;
while n <= numel(t1Ny1)

O1Ny1(n,1) = (t1Ny1(n,1) - tinfNy1(n,1));
O2Ny1(n,1) = (tstNy1 - tinfNy1(n,1));
ONy1(n,1) = (O1Ny1(n,1))/(O2Ny1(n,1));
if (x+.001) > ONy1(n,1) > (x-.001)
    TouNy1(m,1) = tNy1(n,1);
    m = m+1;
end
n = n+1;
end
TouNy1 = TouNy1(1,1);
hNy1 = (pNy1*cNy1*LcNy1)/TouNy1;
BiNy1 = (LcNy1*hNy1)/kNy1;

%% Theoretical Ny1
tinfNy1 = sum(tinfNy1)/numel(tinfNy1);
GrNy1 = (9.81*B*(L^3)*(tstNy1-tinfNy1))/(v^2);
PrNy1 = (c*u)/k;
RaNy1 = PrNy1 * GrNy1;
hconvNy1 = (k/L)*(2+(0.589*RaNy1^(1/4))/(1+(0.469/PrNy1)^(9/16))^(4/9));
hradNy1 = o*e*(tstNy1+tinfNy1)*((tstNy1^2)+(tinfNy1^2));
htotNy1 = hconvNy1+hradNy1;
TouhNy1 = (pNy1*cNy1*LcNy1)/htotNy1;
BihNy1 = (htotNy1*LcNy1)/kNy1;
OhNy1 = exp(-tNy1/TouhNy1);

perdiffTouNy1 = ((TouNy1 - TouhNy1)/TouhNy1)*100; 
perdiffBiNy1 = ((BiNy1 - BihNy1)/BihNy1)*100;

%% Experimental Ny2
n = 1;
m = 1;
while n <= numel(t1Ny2)

O1Ny2(n,1) = (t1Ny2(n,1) - tinfNy2(n,1));
O2Ny2(n,1) = (tstNy2 - tinfNy2(n,1));
ONy2(n,1) = (O1Ny2(n,1))/(O2Ny2(n,1));
if (x+.001) > ONy2(n,1) > (x-.001)
    TouNy2(m,1) = tNy2(n,1);
    m = m+1;
end
n = n+1;
end
TouNy2 = TouNy2(1,1);
hNy2 = (pNy2*cNy2*LcNy2)/TouNy2;
BiNy2 = (LcNy2*hNy2)/kNy2;

%% Theoretical Ny2
tinfNy2 = sum(tinfNy2)/numel(tinfNy2);
GrNy2 = (9.81*B*(D1Ny2^3)*(tstNy2-tinfNy2))/(v^2);
PrNy2 = (c*u)/k;
RaNy2 = PrNy2 * GrNy2;
hconvNy2 = (k/D1Ny2)*(2+(0.589*RaNy2^(1/4))/(1+(0.469/PrNy2)^(9/16))^(4/9));
hradNy2 = o*e*(tstNy2+tinfNy2)*((tstNy2^2)+(tinfNy2^2));
htotNy2 = hconvNy2+hradNy2;
TouhNy2 = (pNy2*cNy2*LcNy2)/htotNy2;
BihNy2 = (htotNy2*LcNy2)/kNy2;
OhNy2 = exp(-tNy2/TouhNy2);

perdiffTouNy2 = ((TouNy2 - TouhNy2)/TouhNy2)*100; 
perdiffBiNy2 = ((BiNy2 - BihNy2)/BihNy2)*100;

%% Experimental Ny3
n = 1;
m = 1;
while n <= numel(t1Ny3)

O1Ny3(n,1) = (t1Ny3(n,1) - tinfNy3(n,1));
O2Ny3(n,1) = (tstNy3 - tinfNy3(n,1));
ONy3(n,1) = (O1Ny3(n,1))/(O2Ny3(n,1));
if (x+.001) > ONy3(n,1) > (x-.001)
    TouNy3(m,1) = tNy3(n,1);
    m = m+1;
end
n = n+1;
end
TouNy3 = TouNy3(1,1);
hNy3 = (pNy3*cNy3*LcNy3)/TouNy3;
BiNy3 = (LcNy3*hNy3)/kNy3;

%% Theoretical Ny3
tinfNy3 = sum(tinfNy3)/numel(tinfNy3);
GrNy3 = (9.81*B*(D1Ny3^3)*(tstNy3-tinfNy3))/(v^2);
PrNy3 = (c*u)/k;
RaNy3 = PrNy3 * GrNy3;
hconvNy3 = (k/D1Ny3)*(2+(0.589*RaNy3^(1/4))/(1+(0.469/PrNy3)^(9/16))^(4/9));
hradNy3 = o*e*(tstNy3+tinfNy3)*((tstNy3^2)+(tinfNy3^2));
htotNy3 = hconvNy3+hradNy3;
TouhNy3 = (pNy3*cNy3*LcNy3)/htotNy3;
BihNy3 = (htotNy3*LcNy3)/kNy3;
OhNy3 = exp(-tNy3/TouhNy3);

perdiffTouNy3 = ((TouNy3 - TouhNy3)/TouhNy3)*100; 
perdiffBiNy3 = ((BiNy3 - BihNy3)/BihNy3)*100;

%% Experimental Steel
n = 1;
m = 1;
while n <= numel(t1Steel)

O1Steel(n,1) = (t1Steel(n,1) - tinfSteel(n,1));
O2Steel(n,1) = (tstSteel - tinfSteel(n,1));
OSteel(n,1) = (O1Steel(n,1))/(O2Steel(n,1));
if (x+.001) > OSteel(n,1) > (x-.001)
    TouSteel(m,1) = tSteel(n,1);
    m = m+1;
end
n = n+1;
end
TouSteel = TouSteel(1,1);
hSteel = (pSteel*cSteel*Lc)/TouSteel;
BiSteel = (Lc*hSteel)/kSteel;

%% Theoretical Steel
tinfSteel = sum(tinfSteel)/numel(tinfSteel);
GrSteel = (9.81*B*(L^3)*(tstSteel-tinfSteel))/(v^2);
PrSteel = (c*u)/k;
RaSteel = PrSteel * GrSteel;
hconvSteel = (k/L)*(2+(0.589*RaSteel^(1/4))/(1+(0.469/PrSteel)^(9/16))^(4/9));
hradSteel = o*0.5*(tstSteel+tinfSteel)*((tstSteel^2)+(tinfSteel^2));
htotSteel = hconvSteel+hradSteel;
TouhSteel = (pSteel*cSteel*Lc)/htotSteel;
BihSteel = (htotSteel*Lc)/kSteel;
OhSteel = exp(-tSteel/TouhSteel);

perdiffTouSteel = ((TouSteel - TouhSteel)/TouhSteel)*100; 
perdiffBiSteel = ((BiSteel - BihSteel)/BihSteel)*100;

%% Making a Table
mat = ['Nylo1';'Nylo2';'Nylo3';'Alumi';'Steel'];
perdiffTou = [perdiffTouNy1; perdiffTouNy2; perdiffTouNy3; perdiffTouAlum; perdiffTouSteel];
perdiffBi = [perdiffBiNy1; perdiffBiNy2; perdiffBiNy3; perdiffBiAlum; perdiffBiSteel];
Bi = [BiNy1; BiNy2; BiNy3; BiAlum; BiSteel];
Bih = [BihNy1; BihNy2; BihNy3; BihAlum; BihSteel];
Tou = [TouNy1; TouNy2; TouNy3; TouAlum; TouSteel];
Touh = [TouhNy1; TouhNy2; TouhNy3; TouhAlum; TouhSteel];
Table = table(mat, Bi, Bih, perdiffBi, Tou, Touh, perdiffTou)

%% Plotting
figure(1);clf
plot(tNy3, ONy3, tNy3, OhNy3)
xlabel('Time (s)')
ylabel('Bi Number')
legend('Measured','Theoretical')

figure(2);clf
plot(tAlum, OAlum, tAlum, OhAlum)
xlabel('Time (s)')
ylabel('Bi Number')
legend('Measured','Theoretical')

figure(3);clf
plot(tNy1, ONy1, tAlum, OAlum, tSteel, OSteel)
xlabel('Time (s)')
ylabel('Bi Number')
legend('Nylon','Aluminium','Steel')

figure(4);clf
plot(tNy1, ONy1, tNy2, ONy2, tNy3, ONy3)
xlabel('Time (s)')
ylabel('Bi Number')
legend('0.5 in','o.75 in','1.0 in')
