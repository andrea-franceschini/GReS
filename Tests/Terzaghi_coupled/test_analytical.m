clear; clc;

%Proprietà materiale (sabbia)
k = 1.0e-5; %[m/s]
porosity = 0.375; 
beta = 4.4e-7 ; % [kPa^-1] fluid compressibility
lambda = 40000; %[kPa] costante lamè
mu = 40000; %[kPa] costante lamè
biot = 1.0; %coefficiente di Biot
E = 100000; %[kPa]
gamma = 9.81; %[kPa] water speicific weight 

%vettore profondità
L = 10; %[m] lunghezza colonna di terreno
z = [0:0.5:L];
z = z';
nz = length(z);
M = (porosity*beta)^-1; %Biot Modulus, assuming cbr=0
Ku = lambda + 2*(mu/3) + biot^2*M;
cm = (lambda+2*mu)^-1; %vertical uniaxial compressibility
c = k/(gamma*(M^-1+biot^2*cm)); %consolidation coefficient

pL = 10; %[kPa] carico esterno 
tau = L^2/c; %[s] charateristic consolidation time
p0 = zeros(length(z),1);
p0(1:end) = (biot*M*pL)/(Ku+4*mu/3);
u0 = arrayfun(@(z) 1/(Ku+4*mu/3)*pL*(L-z),z);

%number of series terms
nm = 10000;
m = [0:nm];

%array for time instants
t = [30 60 600 1200 1800 3600];
nt = length(t);

cellsP = arrayfun(@(m) bsxfun(@(z,t) (1/(2*m+1))*exp((-(2*m+1)^2*pi^2*c*t)/(4*L^2))*sin(((2*m+1)*pi*z)/(2*L)),z,t),m,'UniformOutput',false);
cellsU = arrayfun(@(m) bsxfun(@(z,t) (1/(2*m+1)^2)*exp((-(2*m+1)^2*pi^2*c*t)/(4*L^2))*cos(((2*m+1)*pi*z)/(2*L)),z,t),m,'UniformOutput',false);

matP = cell2mat(cellsP);
matU = cell2mat(cellsU);

seriesP = sum(reshape(matP,nz,nt,[]),3);
seriesU = sum(reshape(matU, nz,nt,[]),3);

p = (4/pi)*seriesP.*p0; %[kPa]
u = cm*((L-z)-(8*L/pi^2)*seriesU).*p0 + u0;

coord = load('coordZ.txt');
ind = zeros(length(coord),1);
for i=1:length(coord)
    ind(i)=find(z==coord(i));
end
u0 = flip(u0);
pIni = p0(ind);
uIni = u0(ind);
save 'ind.dat' ind -ascii 
save 'uIni.dat' uIni -ascii
save 'pIni.dat' pIni -ascii
figure(1)
plot(p,flip(repmat(z,[1,nt])),'-o')


figure(2)
plot(u,flip(repmat(z,[1,nt])),'-o')

%arrayfun(@(m) bsxfun(@(a,b) m+a*b,a,b),m,'UniformOutput',false)


% ptmp(k,j,i) = (4*p0(k))/pi*(1/(2*m(i)+1))*exp((-(2*m(i)+1)^2*pi^2*c*t(j))/(4*L^2))*sin(((2*m(i)+1)*pi*z(k))/(2*L));
%             utmp(k,j,i) = cm*p0(k)*((L-z(k))-(8*L/pi^2)*(1/(2*m(i)+1)^2)*exp((-(2*m(i)+1)^2*pi^2*c*t(j))/(4*L^2))*cos(((2*m(i)+1)*pi*z(k))/(2*L)))+u0(k);
