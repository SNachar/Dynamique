function F = Force(Nn,Nt,dt,type)
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
    case 'creneau';
        F(Nn,:)=type.F0*ones(1,Nt);
    case 'echelon';
        F(Nn,1:type.Ntfin)=type('F0')*ones(1,type.Ntfin);
end
end