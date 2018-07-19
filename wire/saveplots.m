function saveplots(fname)

fname = sprintf('%s',fname);

set(gca,'TitleFontWeight','normal');
set(findall(gcf,'-property','FontSize'),'FontSize',16)
set(findall(gcf,'-property','FontName'),'FontName','Times New Roman')
%set(findall(gcf,'-property','Interpreter'),'Interpreter','Latex')
set(0,'defaultAxesFontName', 'Times New Roman')
set(0,'defaultTextFontName', 'Times New Roman')
print('-dpng','-r300',fname);
print('-dpdf',strrep(fname,'.png','.pdf'));
fprintf('Wrote %s{.png,.pdf}\n',strrep(fname,'.png','.pdf'));