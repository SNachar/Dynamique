function F = Force(Nn,Nt,dt,T,type)
% Fonction Force
% Création d'une matrice contenant en colonne pour chaque pas de temps le
% second membre associé à la force. La force est exercée en bout de poutre.
%
% Nn : Nombre de noeuds
% Nt : Nombre de pas de temps
% dt : Incrément de temps
% type : Obj de type Struct contenant l'information sur le type de
% chargement et l'intensité

F = zeros(Nn,Nt);
switch type.type
    case 'echelon';
        F(Nn,:)=type.F0*ones(1,Nt);
    case 'creneau';
        Ntfin=floor(type.Tfin/dt);
        F(Nn,1:Ntfin)=type.F0*ones(1,Ntfin);
    case 'rampe'
        F (Nn,1:Nt) = (type.F0/ceil(T/dt))*(1:Nt);
    case 'sinus'
        F (Nn,1:Nt) = type.F0*sin((type.w*dt)*(1:Nt));
    otherwise
        error('Mauvais choix de force')
end
end