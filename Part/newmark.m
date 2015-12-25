function [un,vn,an,Sinteg] = newmark(M,K,F,u,v,a,dt,Sinteg)
% R�solution du probl�me � l'aide d'un sch�ma de Newmark
%   M,K,F : matrices du syst�me M*�+K*u=F
%   u : les valeurs de u � l'instant pr�c�dent
%   v : les valeurs de la d�riv�e 1er en temps de u � l'instant pr�c�dent
%   a : les valeurs de la d�riv�e 2nd en temps de u � l'instant pr�c�dent
%   un : les valeurs de u au nouvel instant
%   vn : les valeurs de la d�riv�e 1er en temps de u au nouvel instant
%   an : les valeurs de la d�riv�e 2nd en temps de u au nouvel instant
%   invS : Inverse de la matrice de masse g�n�ralis�e


    % Param�tres du sch�ma de Newmark
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
            warning('Le type est mal d�fini : les seuls possibles sont Newmark, Euler_Forward et Euler_Backward');
            beta1 = 1/2;
            beta2 = 1;
            gamma = 1;
    end
    % Calculs
    if exist('Sinteg') && ~isfield(Sinteg,'invS')
        Sinteg.invS=inv(M+beta2*dt^2*K);
    end
    % Pr�dicateur
    vp = v+dt*(1-gamma)*u;
    up = u+dt*v+dt^2/2*(1-2*beta1)*a;
    % R�sultat
    an = Sinteg.invS*(F-K*up);
    vn = vp + dt*gamma*an;
    un = up + dt^2*beta2*an;
end