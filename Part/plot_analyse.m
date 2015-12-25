function plot_analyse(beta,gamma)
% Extrait de l'analyse spectrale
f = @(Omega,beta,gamma)(abs(Omega.*sqrt(Omega.^2.*beta.*-4.0+Omega.^2.*gamma+Omega.^2.*(1.0./4.0)+Omega.^2.*gamma.^2-4.0).*2.0-Omega.^2.*beta.*4.0+Omega.^2.*gamma.*2.0+Omega.^2-4.0).*(1.0./4.0))./abs(Omega.^2.*beta+1.0);
Omega = 10.^(-3:0.1:3);
    figure1=figure('Name',['Rayon spectral pour \gamma = ',num2str(gamma),', \beta = ',num2str(beta)]);
    axes1 = axes('Parent',figure1);
    hold(axes1,'on');
    xlim(axes1,[0.001 1000]);
    ylim(axes1,[0 1.1]);
    xlabel({'\Omega = h \omega'});
    ylabel({'\rho'});
    title({['Rayon spectral pour \gamma = ',num2str(gamma),', \beta = ',num2str(beta)]});
    box(axes1,'on');
    set(axes1,'XColor',[0 0 0],'XGrid','on','XScale','log','YColor',[0 0 0],...
    'YMinorTick','on','ZColor',[0 0 0]);
        plot(Omega,f(Omega,beta,gamma),'DisplayName',['\beta = ',num2str(beta)]);
end