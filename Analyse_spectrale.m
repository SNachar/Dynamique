syms beta gamma h k m lambda Omega
% Matrice associée à n+1
Mn = [gamma*h*k,m;...
    m+beta*h^2*k,0];
% Matrice associée à n
M = [-(1-gamma)*h*k,m;...
    m-(1/2-beta)*h^2*k,h*m];
% On passe de \dot{u} à h*\dot{u} 
Mn(:,2)=Mn(:,2)./h;
M(:,2)=M(:,2)./h;
Mn = h/m * Mn;
M = h/m * M;
Mn = simplify(subs(Mn,m,k*h^2/Omega^2));
M = simplify(subs(M,m,k*h^2/Omega^2));
% Equation VP généralisée :
eqn = (M - lambda * Mn);
eqn=det(simplify(eqn));
C = coeffs(eqn,lambda);
C=C./C(3);
C=simplify(C);
display('Les coefficients du polynôme caractéristique'); 
pretty(C)
resu=solve(eqn==0,lambda);

f = symfun(abs(resu(1)),[Omega,beta,gamma]);
f2 = symfun(abs(resu(2)),[Omega,beta,gamma]);
hold on
Omega = 10.^(-3:0.1:3);
for gamma = 0.5:0.1:1
    figure1=figure('Name',['Rayon spectral pour \gamma = ',num2str(gamma)]);
    axes1 = axes('Parent',figure1);
    hold(axes1,'on');
    xlim(axes1,[0.001 1000]);
    ylim(axes1,[0.5 1.1]);
    xlabel({'\Omega = h \omega'});
    ylabel({'\rho'});
    title({['Rayon spectral pour \gamma = ',num2str(gamma)]});
    box(axes1,'on');
    set(axes1,'XColor',[0 0 0],'XGrid','on','XScale','log','YColor',[0 0 0],...
    'YMinorTick','on','ZColor',[0 0 0]);
    hold on
    for beta = [0.4, 0.41, 0.4225, 0.5];
        plot(Omega,f2(Omega,beta,gamma),'DisplayName',['\beta = ',num2str(beta)]);
    end
    hold off
end