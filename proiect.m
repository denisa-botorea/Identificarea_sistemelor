t = Botorea(:,1);
u = Botorea(:,2); %semnalul albastru este de intrare
y = Botorea(:,3); %semnalul rosu este de iesire
plot(t, [u ,y]), shg
%% prima metoda de identificare folosind rezonanta
clc
%indecsi alesi pentru semnalul de iesire
i1 =326;  
i2 = 335; 
%indecsi alesi pentru semnalul de intrare
i3 = 323;  
i4 = 333;   


k = mean(y)/mean(u) %eroarea este influentata de k si dt
%amplificarea la rezonanta
Mr = (y(i1)-y(i2))/(u(i3)-u(i4))/k
%perioada la rezonanta pentru semnalul de iesire
Tr1 = 2*(t(i2)-t(i1))  
%Tr2 = 2*(t(i4)-t(i3)) %pentru semnalul de intare
%pulsatia la rezonanta
wr = 2*pi/Tr1 
zeta =sqrt((Mr-sqrt(Mr^2-1))/2/Mr) %tita trebuie sa fie mai mic ca sqrt(2)/2
wn = wr/sqrt(1-2*zeta^2) %legatura dintre wr si wn este prin tita

num = k*wn^2;
den = [1 2*zeta*wn wn^2];
H = tf(num, den)


%spatiul starilor
A = [0,1;-wn^2,-2*zeta*wn];
B = [0;k*wn^2];
C = [1,0];
D = [0]; %transfer instanteu
sys = ss(A,B,C,D);
dt = t(2)-t(1);
yc = lsim(sys,u,t-t(1),[y(1), (y(2)-y(1))/dt]);  %lsim este simularea pe baza modelului identificat
figure
plot(t,[y,yc]), title('yc calculat din spatiul starilor')

%eroarea medie patratica
J = norm(y-yc)/sqrt(length(y))
%eroarea medie patratica normalizata
eMPN = norm(y-yc)/norm(y-mean(y))*100

%% identificarea de pe nyquist
%figure
%nyquist(num, den), title('Nyquist 1')
%bode(num, den), grid
% k este cfj din diagrama nyquist
clc
k = mean(y)/mean(u)

%pentru semnalul de iesire
i5 = 397; %ymax  
i6 = 404; %ymin
%pentru semnalul de intrare
i7 = 393; %umax
i8 = 401; %umin

dt = t(i7)-t(i5); %t(i5)-t(i7)
Tn = 2*(t(i8)-t(i7)) 
wn1 = 2*pi/Tn
ph = (dt*wn1*180)/pi  %ph=90 => wn1=wn
M1 = (y(i5)-y(i6))/(u(i7)-u(i8))/k
zeta1 = k/M1/2
num1 = k*wn1^2;
den1 = [1 2*zeta1*wn1 wn1^2];
% Hnyquist = tf(num1, den1)
% 
% figure
% nyquist(num1, den1), title('Nyquist')

A1 = [0,1;-wn1^2,-2*zeta1*wn1];
B1 = [0;k*wn1^2];
C = [1,0];
D = [0];
sys = ss(A1,B1,C,D);
dt1 = t(2)-t(1)
ycn = lsim(sys,u,t-t(1),[y(1), (y(2)-y(1))/dt1]);  
figure
plot(t,[y,ycn]), title('y calculat din spatiul starilor')

J1 = norm(y-ycn)/sqrt(length(y))
eMPN1 = norm(y-ycn)/norm(y-mean(y))*100


%% identificarea de pe bode
clc
bode(num1, den1), grid %verific panta de -40db/dec
wr1 = wn1*sqrt(1-2*zeta1^2)

% indici pentru semnalul de iesire
i9 = 924; 
i10 = 928; 
%indici pentru semnalul de intrare
i11 = 921; 
i12 = 925; 

%calculul perioadei
T2 = 2*(t(i10)-t(i9)) %trebuie ca T2=Tn/2 , da bine rezultatul
M2 = (y(i9)-y(i10))/(u(i11)-u(i12))/k
amplificarea = 20*log10(M2) 


%% identificare parametrica
clc
% pasul de achizitie
dt = t(2)-t(1);
%datele de identificare
did = iddata(y,u,dt);

%% modelul calculat cu ARX
Marx = arx(did, [2 1 1])  %nA=2, nB=1(nu avem zero), nd=1

%functia de transfer in z/modelul in discret
Hz_arx = tf(Marx.B, Marx.A, dt, 'variable', 'z^-1')
%functia de transfer in s/modelul in continuu
Hs_arx = d2c(Hz_arx, 'zoh') %are un zero si 2 poli in s
%daca scriem tf(Marx.B, Marx.A) obtin functia de transfer in s
Hs = d2c(Hz_arx) %are aceasi forma ca Hs_arx


%validarea statistică
resid(did, Marx, 'corr', 5) %5=numarul de intercorelatii pe care le-am calculat
%gradul de suprapunere
figure; compare(did, Marx) %validarea modelului se face pe aceleasi date de identificare si validare
%intercorelatia este trecuta
%la autocorelatie ne intreseaza nA+nB=3, nu trece testul

%% modelul obtinut prin MVI (metoda variabilelor instrumentale)
clc
MVI = iv4(did, [2 1 1]) %nA=2, nB=1, nd=1

%validare statistică
figure; resid(did, MVI, 'corr', 5)
%gradul de suprapunere
figure; compare(did, Marx, MVI)
%nu trece intercorelatia
%MVI este mai bun decat ARX , uitandu-ne dupa compare

Hz_iv = tf(MVI.B, MVI.A, dt, 'variable', 'z^-1')
Hs_iv = d2c(Hz_iv, 'zoh')

%% modelul obtinut prin metoda output error
clc
Moe = oe(did, [1 2 1]) %nB=1, nF=2, nd=1

%validarea statistica
resid(did, Moe)
%gradul de suprapunere
figure; compare(did, Marx, MVI, Moe)
%Moe este cel mai bun din compare
%intercorelatia trece, autocorelatia nu trece

Hz_oe = tf(Moe.B, Moe.F, dt, 'variable', 'z^-1')
Hs_oe = d2c(Hz_oe, 'zoh')

%% modelul obtinut prin ARMAX
clc
Marmax = armax(did, [2 1 3 1]) %nA=1, nB=1, nC=3, nd=1

%validare statistică
figure; resid(did, Marmax, 'corr', 5) %trec testele de intercorelatie si autocorelatie (esantioanele sunt in banda de incredere)
%gradul de suprapunere
figure; compare(did, Marx, MVI, Moe, Marmax)

Hz_armax = tf(Marmax.B, Marmax.A, dt, 'variable', 'z^-1')
Hs_armax = d2c(Hz_armax, 'zoh')
%% minimizarea erorii de predictie
%modelul de tip spatiul starilor obtinut cu pem
Mpem = pem(did, 2) %ordinul = 2

%validare statistică
resid(did, Mpem, 'corr', 5)
%gradul de suprapunere
figure; compare(did, Mpem)

%% 
clc
% modelul de tip spatiul starilor obtinut cu n4sid
M4 = n4sid(did, 1:10) %retine 10 valori singulare

%validare statistica
figure; resid(did, M4)
%gradul de suprapunere
figure; compare(did, Mpem, M4)

%% rafinarea modelului Marx si MVI (imbunatatirea lor)
clc
Marx_oe = oe(did, Marx)
MVI_oe = oe(did, MVI)

subplot(211), resid(did, Marx_oe), title('Marx_oe')
subplot(212), resid(MVI_oe, did), title('MVI_oe')
figure; compare(did, Marx, Marx_oe, MVI, MVI_oe)

