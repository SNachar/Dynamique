function [dt,Nt]=choix_temps(choix_dt,T,VP)
% Choix de la strat�gie pour dt
% T : Temps total
% VP : Valeurs propres du probl�me (Pour Geradin-Rixen)

switch choix_dt{1}
    case 'Manuel_dt'
        dt = choix_dt {2};
        Nt = ceil(T/dt);
    case 'Manuel_Nt'
        Nt = choix_dt {2};
        dt = T/Nt;
    case 'Geradin-Rixen'
        VPmax=max(VP);
        TPmin=2*pi*VPmax^(-0.5);
        dt = 0.05*TPmin;
        Nt = ceil(T/dt);
    otherwise
        warning(['Non d�finition de la strat�gie pour dt \n', ...
            'Utilisation de Geradin-Rixen']);
        VPmax=max(VP);
        TPmin=2*pi*VPmax^(-0.5);
        dt = 0.05*TPmin;
        Nt = ceil(T/dt);
end
