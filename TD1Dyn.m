%% TD 1 de Dynamique des Structures
% S.Nachar & M.Tirvaudey (2015)
%
% Etude des schémas d'intégration temporelle via le cas d'une poutre 
% encastrée en traction-compression

%% Mise à zéro :
clearvars;
clc;
close all;
addpath('./Part');

%% Paramètres :

%%% Paramètres variables
Ne = 1000; % Nombre d'éléments
T = 0.005; %Durée totale

%%%% Choix de la stratégie pour l'incrément de temps
% choix_dt = {'Manuel_dt',1E-4}; % Choix de l'incrément de temps
% choix_dt = {'Manuel_Nt',1000}; % Choix du pas de temps
choix_dt = {'Geradin-Rixen'}; % dt = 0.05*TPmin;

%%%% Choix du schéma d'intégration
% Sinteg = struct('type','Euler_Backward');
% Sinteg = struct('type','Euler_Forward');
% Sinteg = struct('type','Newmark','coeff',[0 1/2]);
 Sinteg = struct('type','Newmark_alpha','alpha',0.5); 

%%%% Force :
 force = struct('type','creneau','F0',30,'Tfin',0.5);
% force = struct('type','echelon','F0',30);
% force = struct('type','rampe','F0',30);
% force = struct('type','sinus','F0',30,'w',1E1);

 
%%%% Troncature des modes propres
Redu = true; % Réduction du modèle
Nm = 20; %Nombre de modes gardés


%%% Paramètres de l'énoncé
L = 1; % Longueur totale
Rho = 70E3; % Masse volumique
E = 220E9; % Module de Young
A = 1E-4; % Section
% Coefficients de Kv :
a = 0;
b = 0;

% Paramètres implicites :
Nn = Ne + 1; % Nombre de noeuds

%% Creation des matrices de base :

le = L/Ne;
k = E*A/(2*le);
m = Rho * A * le /6;

Ke=k*[1 -1; -1 1]; % matrice de rigidite
Me=m*[2 1; 1 2]; % et de masse elementaire

%% Assemblage
% Matrice de masse & de rigidité :
K = zeros(Nn,Nn);
M = zeros(Nn,Nn);
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
    VP=sort(nonzeros(D));
% trie des v.propres
if Redu
    [Vs,Id]=sort(diag(D));
% Arrangement de V
    W = V(:,Id);
% Conservation des premiers modes propres
    Q = W(1:end,1:Nm);
    VP=VP(1:Nm);
end
% Définition de l'incrément de temps
    [dt,Nt]=choix_temps(choix_dt,T,VP);

%% Résolution

F = Force(Nn,Nt,dt,T,force);
U = zeros(Nn,Nt);
Up = zeros(Nn,Nt);
Upp = zeros(Nn,Nt);

% Conditions Limites (Substitution)
[M,K,F,U,Up,Upp]=CL(M,K,F,U,Up,Upp);
if Redu; Q(1,:)=[]; end

%Changement de base des matrices
if Redu
    K=Q'*K*Q;    M=Q'*M*Q;    F=Q'*F(:,1:Nt-1);
    U=Q'*U(:,1:Nt-1);    Up=Q'*Up(:,1:Nt-1);    Upp=Q'*Upp(:,1:Nt-1);
end

% Intégration temporelle
for it=1:Nt-1
   [U(:,it+1),Up(:,it+1),Upp(:,it+1),Sinteg]... 
         = newmark(M,K,F(:,it),U(:,it),Up(:,it),Upp(:,it),dt,Sinteg);
end

% Changement de base des vecteurs
if Redu
    U=Q*U;    Up=Q*Up;    Upp=Q*Upp;
end

% Condition limite en x=0 replacée
U=[zeros(1,Nt);U];  Up=[zeros(1,Nt);Up];    Upp=[zeros(1,Nt);Upp];

% Figures
if strcmp(Sinteg.type,'Newmark_alpha')
    gamma = 0.5+Sinteg.alpha;
    beta = 0.25*(gamma+0.5)^2;
    plot_analyse(beta,gamma)
elseif strcmp(Sinteg.type,'Newmark')
    plot_analyse(Sinteg.coeff(1),Sinteg.coeff(2))
end
anim_poutre(U,Nn,Nt,force,Sinteg)