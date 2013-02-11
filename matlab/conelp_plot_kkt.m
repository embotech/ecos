function conelp_plot_kkt(K,name)

figure; subplot(1,2,1); spy(K,'k.'); title(name); 
[LK,DK,PK] = ldl(K); 
subplot(1,2,2); spy(LK,'k.'); title(['L',name]);
savefig(['fillin_',name],gcf,'pdf')