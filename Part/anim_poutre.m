function anim_poutre(U,Nn,Nt,force,Sinteg)

figure1 = figure;
axes1 = axes('Parent',figure1);
y=U(:,1);
h = plot(y,'YDataSource','y');
xlabel('N° Noeud');
ylabel('Ux (m)');
bylim=[1.1*min(U(Nn,:)) 1.1*max(U(Nn,:))];
ylim(axes1,bylim);
id = 0;
for k = 2:ceil(Nt/1000):Nt
   id = id+1;
   y = U(:,k);
   refreshdata(h,'caller')
   ylim(axes1,bylim);
   drawnow
   if mod(id,100)==0
   print(['Images/',force.type,'-',Sinteg.type,...
            '-k=',num2str(k),'.png'],'-dpng');
   end
end