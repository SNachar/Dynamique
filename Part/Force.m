function F = Force(Nn,Nt,dt,T,type)
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