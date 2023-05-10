clc;
close all;
fh = figure();
fh.WindowState = 'maximized';

data=load("Comsol\LineGraphU.txt");
x1 = data(:,1);
y1 = data(:,2);
data = load('C:\Users\gpesc\Desktop\GitProg\KeyReal\main\result.txt');
x2 = data(:,1)/1000;
y2 = data(:,2)/1000;

plot(x1, y1, 'o','markerSize', 15);
hold on;
plot(x2, y2,'LineWidth', 3);
legend("Обратное преобразование");
set(gca, 'FontSize',36);
grid on;
set(gcf,'color','w');
set(gcf,'color','w');

legend('Comsol','Интеграл','Location','northeast');
%title('$\mathbf{u}^+_{ss}$', 'Interpreter', 'Latex', 'FontSize',75);
xlabel('$x$', 'FontSize',16, 'Interpreter', 'LaTeX');
%set(get(gca,'title'),'Position',[0 0])
set(gca, 'GridAlpha', 1);

set(findall(gcf,'type','text'), 'FontSize', 65, 'FontName', 'Times New Roman')
set(gca,'FontSize', 48, 'FontName', 'Times New Roman')
pbaspect([1.62 1 1])
set(gcf,'units','centimeters','position',[0,0,6.2,3.8])