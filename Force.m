function F = Force(Nn,Nt,dt,type)
% Fonction Force
% Cr�ation d'une matrice contenant en colonne pour chaque pas de temps le
% second membre associ� � la force. La force est exerc�e en bout de poutre.
%
% Nn : Nombre de noeuds
% Nt : Nombre de pas de temps
% dt : Incr�ment de temps
% type : Obj de type Struct contenant l'information sur le type de
% chargement et l'intensit�

F = zeros(Nn,Nt);
switch type.type
    case 'creneau';
        F(Nn,:)=type.F0*ones(1,Nt);
    case 'echelon';
        F(Nn,1:type.Ntfin)=type('F0')*ones(1,type.Ntfin);
end
end