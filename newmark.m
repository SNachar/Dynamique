function [u,v,a,Sinteg] = newmark(M,K,F,un,vn,an,dt,Sinteg)
% Résolution du problème à l'aide d'un schéma de Newmark
%   M,K,F : matrices du système M*ü+K*u=F
%   un : les valeurs de u à l'instant précédent
%   vn : les valeurs de la dérivée 1er en temps de u à l'instant précédent
%   an : les valeurs de la dérivée 2nd en temps de u à l'instant précédent
%   u : les valeurs de u au nouvel instant
%   v : les valeurs de la dérivée 1er en temps de u au nouvel instant
%   a : les valeurs de la dérivée 2nd en temps de u au nouvel instant
%   invS : Inverse de la matrice de masse généralisée

    % Paramètres du schéma de Newmark
    if exist('Sinteg')
        beta1 = Sinteg.coeff(1);
        beta2 = Sinteg.coeff(2);
        gamma = Sinteg.coeff(3);
    else
        beta1 = 0;
        beta2 = 1/2;
        gamma = 1;
    end
    
    % Calculs
    if exist('Sinteg') || ~isfield(Sinteg,'invS')
        Sinteg.invS=(M+dt^2*K*(gamma+beta1));
    end
    a = Sinteg.invS\(F-K*(un+dt*vn+dt^2*(3/2-gamma-beta1)*an));
%    v = zeros(size(K,2),1); % TODO calcul de la dérivée 1er
    v = vn + dt*((1-gamma)*an + gamma*a);
%    u = zeros(size(K,2),1); % TODO calcul du champs
    u = un + dt*v + dt^2*((1/2-beta2)*an + beta1*a);
end