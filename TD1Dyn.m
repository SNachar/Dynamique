%% TD 1 Dynamique

%% Paramètres :
Ne = 20; % Nombre d'éléments
T = 0.01; %Durée totale
Nt = 1000; % Incrément de temps
L = 1; % Longueur totale
Rho = 70E3; % Masse volumique
E = 220E9; % Module de Young
A = 1E-4; % Section
Nm = 20; %Nombre de modes gardés
% Coefficients de Kv :
a = 0;
b = 0;
% Force :
type = struct('type','creneau','F0',30);
% Choix du schéma d'intégration
Sinteg = struct('type','Newmark','coeff',[0 1/2 1]);

% Paramètres implicites :
Nn = Ne + 1; % Nombre de noeuds
dt = T/Nt;

%% Creation des matrices de base :

le = L/Ne;
k = E*A/(2*le);
m = Rho * A * le /6;

Ke=k*[1 -1; -1 1]; % matrice de rigidite
Me=m*[2 1; 1 2]; % et de masse elementaire

%% Assemblage
% Matrice de masse & de rigidité :
K = zeros(Nn);
M = zeros(Nn);
Kv = a*K+b*M;
for i=1:Ne
    for j=0:1
        for k=0:1
            K(i+j,i+k)=K(i+j,i+k)+Ke(1+j,1+k);
            M(i+j,i+k)=M(i+j,i+k)+Me(1+j,1+k);
        end;
    end;
end;

%% Recherche de modes propres
% resolution mode propre
[V,D]=eig(M\K);
% trie des v.propres
[Vs,Id]=sort(diag(D));
% Arrangement de V
W = V(:,Id);

% Conservation des premiers modes propres
Q=zeros(Ne+1,Nm);
for i=1:Nm
    Q(:,i)=W(:,i);
end
%changement des valeurs des matrices
% Kq=Q'*K*Q;
% Mq=Q'*M*Q;
% Uq=Q'*U ;
% Fq=Q'*Fcren ;


%% Résolution Euler

F = Force(Nn,Nt,dt,type);
U = zeros (Nn,Nt);
Up = zeros (Nn,Nt);
Upp = zeros (Nn,Nt);
[M,K,F,U,Up,Upp]=CL(M,K,F,U,Up,Upp);
for it=1:Nt-1
   [U(:,it+1),Up(:,it+1),Upp(:,it+1),Sinteg]... 
         = newmark(M,K,F(:,it),U(:,it),Up(:,it),Upp(:,it),dt,Sinteg);
end

U=[zeros(1,Nt);U];
Up=[zeros(1,Nt);Up];
Upp=[zeros(1,Nt);Upp];

figure1 = figure;
axes1 = axes('Parent',figure1);
y=U(:,1)+dt*[0:Nn-1]';
h = plot(y,'YDataSource','y');
for k = 2:Nt
   y = U(:,k);
   ylim(axes1,[0 1.1*U(Nn,k)]);
   ylabel=strcat('k=',num2str(k));
   refreshdata(h,'caller')
   drawnow
end