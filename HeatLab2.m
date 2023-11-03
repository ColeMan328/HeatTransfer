% Cole Stoumbaugh
clear;clc

%% Variable Setup
Test1 = [3.000000 0.050000 0.780000 0.005000 20.814154 0.252798 42.860338 0.255796 36.597565 0.242975 31.771338 0.256966 28.836944 0.264190 27.022422 0.255403 25.580637 0.262274];
Test2 = [4.000000 0.050000 1.060000 0.005000 20.895776 0.276069 57.317171 0.260765 46.595844 0.262932 38.455297 0.269783 33.545760 0.277061 30.314683 0.257371 27.860241 0.262137];
Test3 = [5.000000 0.050000 1.340000 0.005000 20.878864 0.280474 76.695123 0.259732 60.019266 0.247171 47.386184 0.271362 39.815271 0.275409 34.847832 0.258001 31.147913 0.260403];

v = [Test1(1,1), Test2(1,1), Test3(1,1)];
dv = [Test1(1,2), Test2(1,2), Test3(1,2)];
i = [Test1(1,3), Test2(1,3), Test3(1,3)];
di = [Test1(1,4), Test2(1,4), Test3(1,4)];
Tinf = [Test1(1,5), Test2(1,5), Test3(1,5)] + 273;
dTinf = [Test1(1,6), Test2(1,6), Test3(1,6)];
T1 = [Test1(1,7), Test2(1,7), Test3(1,7)] + 273;
dT1 = [Test1(1,8), Test2(1,8), Test3(1,8)];
T2 = [Test1(1,9), Test2(1,9), Test3(1,9)] + 273;
dT2 = [Test1(1,10), Test2(1,10), Test3(1,10)];
T3 = [Test1(1,11), Test2(1,11), Test3(1,11)] + 273;
dT3 = [Test1(1,12), Test2(1,12), Test3(1,12)];
T4 = [Test1(1,13), Test2(1,13), Test3(1,13)] + 273;
dT4 = [Test1(1,14), Test2(1,14), Test3(1,14)];
T5 = [Test1(1,15), Test2(1,15), Test3(1,15)] + 273;
dT5 = [Test1(1,16), Test2(1,16), Test3(1,16)];
T6 = [Test1(1,17), Test2(1,17), Test3(1,17)] + 273;
dT6 = [Test1(1,18), Test2(1,18), Test3(1,18)];
T = [T1;T2;T3;T4;T5;T6];
dT = [dT1;dT2;dT3;dT4;dT5;dT6];
Ts = [(sum(T(1:6,1))/6),(sum(T(1:6,2))/6),(sum(T(1:6,3))/6)]; 
dTs = [(sqrt(sum(dT(1:6,1).^2))/6),(sqrt(sum(dT(2:6,1).^2))/6),(sqrt(sum(dT(3:6,1).^2))/6)]; 

D = 0.375*.0254;
L = 14*.0254;
k = 121;
e = 1;
o = 5.67*10^(-8);
s = 2;

x = linspace(0,L,1000);

%% Calculations

m = 1;
for m = 1:3
    h(1,m) = 1.32*((Ts(1,m)-Tinf(1,m))/D)^0.25 + (o*e*(Ts(1,m)+Tinf(1,m))*(Ts(1,m)^2 + Tinf(1,m)^2));   % Heat Trans Coeff
    dh1(1,m) = ((1.32/(D^.25))*(.25/(Ts(1,m)-Tinf(1,m))^.75) + o*e*((3*Ts(1,m)^2)+(Tinf(1,m)^2)+(2*Ts(1,m)*Tinf(1,m))))*dTs(1,m);   % Heat Trans Coeff
    dh2(1,m) = ((1.32/(D^.25))*(-.25/(Ts(1,m)-Tinf(1,m))^.75) + o*e*((3*Tinf(1,m)^2)+(Ts(1,m)^2)+(2*Ts(1,m)*Tinf(1,m))))*dTinf(1,m);
    dh(1,m) = sqrt(dh1(1,m)^2 + dh2(1,m)^2);
    M(1,m) = sqrt((4*h(1,m))/(k*D));    % m calc
    dM(1,m) = (h(1,m)^(1/2)) * sqrt(1/(k*D)) * dh(1,m); % m calc uncert
    mL(1,m) = M(1,m)*L; % mL calc
    qb(1,m) = k*(pi/4)*(D^2)*(T(1,m)-Tinf(1,m)) * (M(1,m)*tanh(mL(1,m)));   % heat transf at base of fin
    dq1(1,m) = (-k*(pi/4)*(D^2)*(Tinf(1,m)) * (M(1,m)*tanh(mL(1,m)))) * dT(1,m);    % heat transf at base of fin uncert
    dq2(1,m) = (k*(pi/4)*(D^2)*(T(1,m)) * (M(1,m)*tanh(mL(1,m)))) * dTinf(1,m);
    dq3(1,m) = (-k*(pi/4)*(D^2)*(T(1,m)-Tinf(1,m)) * (-tanh(mL(1,m)) - mL(1,m)*sech(mL(1,m)^2))) * dM(1,m);
    dqb(1,m) = sqrt((dq1(1,m)^2)+(dq2(1,m)^2)+(dq3(1,m)^2));
    qtot(1,m) = v(1,m)*i(1,m);  % total heat transf
    dqtot(1,m) = sqrt((dv(1,m)*i(1,m))^2 + (v(1,m)*di(1,m))^2); % total heat transf uncert
    pDiff(1,m) = (qb(1,m)-qtot(1,m))/qtot(1,m);
end

%% Calculating Temp

Tx1 = Tinf(1,1) + ((T(1,1)-Tinf(1,1))*(cosh(M(1,1)*(L-x))/cosh(mL(1,1))));
Tx2 = Tinf(1,2) + ((T(1,2)-Tinf(1,2))*(cosh(M(1,2)*(L-x))/cosh(mL(1,2))));
Tx3 = Tinf(1,3) + ((T(1,3)-Tinf(1,3))*(cosh(M(1,3)*(L-x))/cosh(mL(1,3))));

%% Making Plot

figure(1);clf
plot(x,Tx1,'b');hold on
plot(x,Tx2,'r');hold on
plot(x,Tx3,'g')
legend('Case 1', 'Case 2', 'Case 3')
xlabel('Position')
ylabel('Temperature')

%% Making Table

Case = ['Test 1';'Test 2';'Test 3'];
h_tot = [h(1,1);h(1,2);h(1,3)];
dh_tot = [dh(1,1);dh(1,2);dh(1,3)];
q_base = [qb(1,1);qb(1,2);qb(1,3)];
dq_base = [dqb(1,1);dqb(1,2);dqb(1,3)];
q_tot = [qtot(1,1);qtot(1,2);qtot(1,3)];
dq_tot = [dqtot(1,1);dqtot(1,2);dqtot(1,3)];
percentDiff = [pDiff(1,1);pDiff(1,2);pDiff(1,3)];
tab = table(Case, h_tot, dh_tot, q_base, dq_base, q_tot, dq_tot, percentDiff)
