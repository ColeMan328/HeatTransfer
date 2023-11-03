%Cole Stoumbaugh
clear;clc

%% Setting Variables
%Lablels = ['V [V]',	'DeltaV [V]',	'I [A]',	'DeltaI [A]',	'T1 [deg C]',	'DeltaT1 [deg C]',	'T2 [deg C]',	'DeltaT2 [deg C]',	'T3 [deg C]',	'DeltaT3 [deg C]',	'T4 [deg C]',	'DeltaT4 [deg C]',	'T5 [deg C]',	'DeltaT5 [deg C]',	'T6 [deg C]',	'DeltaT6 [deg C]'];
b1 = [9.82 0.005 1.05 0.005 83.129186 0.26709 80.672668 0.250944 78.572523 0.261683 69.239409 0.256651 62.16139 0.31289 54.8923 0.24637];
b2 = [9.97 0.005 1.03 0.005 55.685855 0.276447 52.501237 0.278431 50.33253 0.340071 44.348711 0.278118 41.983859 0.272884 40.35306 0.285346];
a = [9.84 0.005 1.02 0.005 47.03841 0.280013 44.096006 0.306885 41.652405 0.276545 38.521994 0.279011 37.210754 0.292744 35.782109 0.271199];
l = [.015,.015,.015];
x3 = linspace(0,.001,.075);
x = [0,.015,.030,.045,.06,.075];
x1 = [0,.015,.030, .0375];
x2 = [.0375, .045,.06,.075];
k1 = 121;
k = [121,180,121];
d1 = .025;
d2 = [.025,.025,.013];
v = [b1(1,1), b2(1,1), a(1,1)];
dv = [b1(1,2), b2(1,2), a(1,2)];
i = [b1(1,3), b2(1,3), a(1,3)];
di = [b1(1,4), b2(1,4), a(1,4)];
t1 = [b1(1,5), b2(1,5), a(1,5)];
dt1 = [b1(1,6), b2(1,6), a(1,6)];
t2 = [b1(1,7), b2(1,7), a(1,7)];
dt2 = [b1(1,8), b2(1,8), a(1,8)];
t3 = [b1(1,9), b2(1,9), a(1,9)];
dt3 = [b1(1,10), b2(1,10), a(1,10)];
t4 = [b1(1,11), b2(1,11), a(1,11)];
dt4 = [b1(1,12), b2(1,12), a(1,12)];
t5 = [b1(1,13), b2(1,13), a(1,13)];
dt5 = [b1(1,14), b2(1,14), a(1,14)];
t6 = [b1(1,15), b2(1,15), a(1,15)];
dt6 = [b1(1,16), b2(1,16), a(1,16)];
Tb1 = [t1(1,1),t2(1,1),t3(1,1),t4(1,1),t5(1,1),t6(1,1)];
Tb2 = [t1(1,2),t2(1,2),t3(1,2),t4(1,2),t5(1,2),t6(1,2)];
Ta = [t1(1,3),t2(1,3),t3(1,3),t4(1,3),t5(1,3),t6(1,3)];
dTb1 = [dt1(1,1),dt2(1,1),dt3(1,1),dt4(1,1),dt5(1,1),dt6(1,1)];
dTb2 = [dt1(1,2),dt2(1,2),dt3(1,2),dt4(1,2),dt5(1,2),dt6(1,2)];
dTa = [dt1(1,3),dt2(1,3),dt3(1,3),dt4(1,3),dt5(1,3),dt6(1,3)];

%% Calculations
n = 1;
% best fit line
for n = 1:3
    m = 1;
for m = 1:4
T1(n,m) = ((t3(1,n)-t2(1,n))/l(1,n))*x1(1,m) + 2*t2(1,n) - t3(1,n);
T2(n,m) = ((t5(1,n)-t4(1,n))/l(1,n))*x2(1,m) + t4(1,n) - 3*(t5(1,n)-t4(1,n));
m = m+1;
end
n = n+1;
end

% border temp
TB1 = (3/2).*t3 - .5.*t2;
TB2 = (3/2).*t4 - .5.*t5;
Tdiff = (3/2).*t3 - .5.*t2 - (3/2).*t4 + .5.*t5;
% heat rate
q = v.*i;
dq = sqrt((v.*di).^2 + (i.*dv).^2);
qr = -1 .* k1 .* pi .* (d1^2) .* (t3-t2) .* (.25 ./ l);
dqr = sqrt((.25 .* k .* pi .* (d1.^2) .* (1./l).*dt2).^2 + (.25 .* k .* pi .* (d1.^2) .* (1./l).*dt3).^2);
% Thermal conductivity
k2 = -(4 .* qr .* l)./(pi .* (d2.^2) .* (t5 - t4));
%dk2 = sqrt(((4 .* dqr .* l)/(pi .* d2.^2 .* (t5-t4))).^2 + ((4 .* qr .* t4 .* l)./(pi .* d2.^2 .* (t5-t4).^2 .* t5)).^2 + ((4 .* qr .* l)/(pi .* d2.^2 .* (t5-t4).^2 .* t4)).^2);
dk3 = k2 .* sqrt((dqr./qr).^2 + (dt4./(t5-t4)).^2 + (dt5./(t5-t4)).^2);
%% Plotting
figure(1);clf
errorbar(x, Tb1, dTb1,'.b');hold on
errorbar(x, Tb2, dTb2,'.r');hold on
errorbar(x, Ta, dTa,'.g');hold on
plot(x1,T1(1,:),'b');hold on
plot(x1,T1(2,:),'r');hold on
plot(x1,T1(3,:),'g');hold on
plot(x2,T2(1,:),'b');hold on
plot(x2,T2(2,:),'r');hold on
xlabel('Position (m)');
plot(x2,T2(3,:),'g');
ylabel('Temperature (C)');
%title('Rod Temperature');
legend('Brass 1','Brasss 2','Alumin');


%% Making Chart
Section = ['Brass 1'; 'Brass 2'; "Aluminium"];
Temp1 = reshape(TB1,[3,1]);
Temp2 = reshape(TB2,[3,1]);
TempDiff = reshape(Tdiff,[3,1]);
ThermalCond = reshape(k2, [3,1]);
ThermalCondUncert = reshape(dk3, [3,1]);
HeatTransElec = reshape(q, [3,1]);
HeatTransElecUncert = reshape(dq, [3,1]);
HeatTransTemp = reshape(qr, [3,1]);
HeatTransTempUncert = reshape(dqr, [3,1]);
BoarderTemp = table(Section, Temp1, Temp2, TempDiff, ThermalCond, ThermalCondUncert, HeatTransElec, HeatTransElecUncert, HeatTransTemp, HeatTransTempUncert)
