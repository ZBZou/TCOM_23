%% savefig
function mysavefig(filename)

name1 = append(filename, '.fig');
name2 = append(filename, '.pdf');
saveas(gca, name1);
exportgraphics(gca, name2);
