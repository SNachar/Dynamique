function [un,vn,an,Sinteg] = newmark(M,K,F,u,v,a,dt,Sinteg)
% Résolution du problème à l'aide d'un schéma de Newmark
%   M,K,F : matrices du système M*ü+K*u=F
%   u : les valeurs de u à l'instant précédent
%   v : les valeurs de la dérivée 1er en temps de u à l'instant précédent
%   a : les valeurs de la dérivée 2nd en temps de u à l'instant précédent
%   un : les valeurs de u au nouvel instant
%   vn : les valeurs de la dérivée 1er en temps de u au nouvel instant
%   an : les valeurs de la dérivée 2nd en temps de u au nouvel instant
%   invS : Inverse de la matrice de masse généralisée


    % Paramètres du schéma de Newmark
    switch Sinteg.type
        case 'Newmark'
            beta1 = Sinteg.coeff(1);
            beta2 = Sinteg.coeff(1);
            gamma = Sinteg.coeff(2);
        case 'Newmark_alpha'
            gamma = 0.5+Sinteg.alpha;
            beta1 = 0.25*(gamma+0.5)^2;
            beta2 = beta1;
        case 'Euler_Backward'
            beta1 = 1/2;
            beta2 = 1;
            gamma = 1;
        case 'Euler_Forward'
            beta1 = 1/2;
            beta2 = 0;
            gamma = 0;
        otherwise
            warning('Le type est mal défini : les seuls possibles sont Newmark, Euler_Forward et Euler_Backward');
            beta1 = 1/2;
            beta2 = 1;
            gamma = 1;
    end
    % Calculs
    if exist('Sinteg') && ~isfield(Sinteg,'invS')
        Sinteg.invS=inv(M+beta2*dt^2*K);
    end
    % Prédicateur
    vp = v+dt*(1-gamma)*u;
    up = u+dt*v+dt^2/2*(1-2*beta1)*a;
    % Résultat
    an = Sinteg.invS*(F-K*up);
    vn = vp + dt*gamma*an;
    un = up + dt^2*beta2*an;
end