t = Botorea(:,1);
u = Botorea(:,2); %semnalul albastru este de intrare
y = Botorea(:,4); %semnalul rosu este de iesire
plot(t, [u ,y]), shg

%% identificarea prin fenomenul de rezonanta
clc
%indecsi alesi pentru semnalul de iesire
i1 = 379; %ymax  
i2 = 387; %ymin  
%indecsi alesi pentru semnalul de intrare
i3 = 377; %umax  
i4 = 385; %umin  

%factorul de proporționalitate în regim staționar
k = mean(y)/mean(u)
%amplificarea la rezonanță
Mr = (y(i1)-y(i2))/(u(i3)-u(i4))/k
%perioada la rezonanță pentru semnalul de intrare
Tr = 2*(t(i4)-t(i3));
%pulsația la rezonanță
wr = 2*pi/Tr
zeta =sqrt((Mr-sqrt(Mr^2-1))/2/Mr)
wn = wr/sqrt(1-2*zeta^2)


%faza la rezonanță
phr = (t(i3)-t(i1))*wr; %da in radiani
phr = (t(i3)-t(i1))*wr*180/pi; %da in grade
%constanta de timp a zeroului 
Tz = tan(phr + atan(sqrt(1 - 2*zeta^2)/zeta))/wr

num = k*wn^2*[Tz, 1];
den = [1, 2*zeta*wn, wn^2];
H = tf(num, den)

%spațiul stărilor
A = [0, 1; -wn^2, -2*zeta*wn];
g1 = k*Tz*wn^2;
g2 = k*wn^2-k*Tz*wn^3*2*zeta;
B = [g1 ;g2];
C = [1 0];
D = [0];
dt = t(2)-t(1);
yc = lsim(A,B,C,D, u,t-t(1), [y(1), (y(2)-y(1))/dt - g1*u(1)]);
plot(t, [y yc]), title('y calculat din spatiul starilor')

%eroarea medie patratica
J = norm(y-yc)/sqrt(length(y))
%eroarea medie patratica normalizata
eMPN = norm(y-yc)/norm(y-mean(y))*100




%% identificarea parametrica pentru alte date de identificare 
idx = [78:665];
dt = t(2)-t(1);
did = iddata(y(idx,1), u(idx,1), dt);


%% modelul calculat cu ARX pentru did
clc
Marx = arx(did, [2 2 1])  %nA=2, nB=2(avem un zero), nd=1

%validarea statistică
resid(did, Marx) 
%gradul de suprapunere
figure; compare(did, Marx) 

%% modelul obtinut prin MVI pentru did
clc
MVI = iv4(did, [2 2 1]) %nA=2, nB=2, nd=1

%validare statistică
figure; resid(did, MVI) %trece intercorelatia
%gradul de suprapunere
figure; compare(did, Marx, MVI)

Hz_iv = tf(MVI.B, MVI.A, dt, 'variable', 'z^-1')
Hs_iv = d2c(Hz_iv, 'zoh')

%% modelul obtinut prin metoda output error pentru did
Moe = oe(did, [2 2 1]) %nB=2, nF=2, nd=1

%validarea statistica
resid(did, Moe) %nu trece niciun test
%gradul de suprapunere
figure; compare(did, Marx, MVI, Moe)
Hz_oe = tf(Moe.B, Moe.F, dt, 'variable', 'z^-1')
Hs_oe = d2c(Hz_oe, 'zoh')



%% modelul obtinut prin ARMAX  pentru did
Marmax = armax(did, [2 2 5 1]) %nA=1, nB=2, nC=5, nd=1

%validare statistică
figure; resid(did, Marmax, 'corr', 5) %trec amandoua
%gradul de suprapunere
figure; compare(did, Marx, MVI, Moe, Marmax)

%% modelul obtinut cu pem di n4sid pentru did
Mpem = pem(did, 2) %ordinul = 2
%modelul obtinut cu n4sid
M4 = n4sid(did, 1:15) %retine 15 valori singulare

%validare statistică
subplot(211), resid(did, Mpem, 'corr', 5), title('Mpem') %trece doar intercorelatia
subplot(212), resid(did, M4, 'corr', 5), title('M4') %nu trece nimic
%gradul de suprapunere
figure; compare(did, Marx, MVI, Moe, Marmax, Mpem, M4)

%% rafinarea modelului Marx si M4(imbunatatirea lui )
Marx_oe = oe(did, Marx)
M4_pem = pem(did, M4)

subplot(311), resid(did, Marx_oe, 'corr', 5), title('Marxoe')
subplot(312), resid(did, M4_pem, 'corr', 5), title('M4pem')
figure; compare(did, Marx, Marx_oe, M4, M4_pem)


