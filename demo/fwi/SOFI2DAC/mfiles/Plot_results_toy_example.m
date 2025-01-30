% This is a matlab script to plot the results of the toy example of the DENISE code.

clear all
close all
clc

iteration_no=137;

% general informations about the model
nx=252;		% grid points in x direction
ny=75;		% grid points in y direction
DH=0.2;		% spatial sampling interval
X=(1:1:nx)*DH;	% definition of x axis for plots
Y=(1:1:ny)*DH;	% definition of y axis for plots
FW=20;		% width of PML frame in gridpoints
sources=[5.0, 15.0, 25.0, 35.0, 45.0];	% x coordinate of source positions

% general settings for plots information
fontsize=14;
x_line=[FW FW nx-FW nx-FW]*DH;
y_line=[1 ny-FW ny-FW 1]*DH;

% reading the true model
true_vs=read_binary_matrix(nx,ny,['../par/model/mod_toy_example_true_vs_it0.bin']);
true_vp=read_binary_matrix(nx,ny,['../par/model/mod_toy_example_true_vp_it0.bin']);
true_rho=read_binary_matrix(nx,ny,['../par/model/mod_toy_example_true_rho_it0.bin']);

% reading initial model
start_vs=read_binary_matrix(nx,ny,['../par/model/toy_example/mod_toy_example_vs_it0.bin']);
start_vp=read_binary_matrix(nx,ny,['../par/model/toy_example/mod_toy_example_vp_it0.bin']);
start_rho=read_binary_matrix(nx,ny,['../par/model/toy_example/mod_toy_example_rho_it0.bin']);

% reading inverted vs model
inv_vs=read_binary_matrix(nx,ny,['../par/model/toy_example/mod_toy_example_vs_it' int2str(iteration_no) '.bin']);

% finding minimum and maximum values for colorbar
max_vs=max([max(true_vs(:)) max(start_vs(:)) max(inv_vs(:))]);
min_vs=min([min(true_vs(:)) min(start_vs(:)) min(inv_vs(:))]);

max_vp=max([max(true_vp(:)) max(start_vp(:))]);
min_vp=min([min(true_vp(:)) min(start_vp(:))]);

max_rho=max([max(true_rho(:)) max(start_rho(:))]);
min_rho=min([min(true_rho(:)) min(start_rho(:))]);

% Plotting true model
figure('units','normalized','outerposition',[0.75 0 0.25 1])
subplot(3,2,1)
imagesc(X,Y,true_vs)
hold on
line(x_line,y_line,'LineStyle','--','Color','k')
plot(sources,0.2*ones(1,5),'*r','MarkerSize',12)
hold off
xlabel('x in m','FontSize',fontsize)
ylabel('y in m','FontSize',fontsize)
title(['true S-wave velocity model in m/s'],'FontSize',fontsize)
colb=colorbar;
coll=get(colb,'ylabel');
set(coll,'String','vs in m/s'); 
caxis([min_vs max_vs])
axis equal tight

subplot(3,2,2)
imagesc(X,Y,true_vp)
hold on
line(x_line,y_line,'LineStyle','--','Color','k')
plot(sources,0.2*ones(1,5),'*r','MarkerSize',12)
hold off
xlabel('x in m','FontSize',fontsize)
ylabel('y in m','FontSize',fontsize)
title(['true P-wave velocity model in m/s'],'FontSize',fontsize)
colb=colorbar;
coll=get(colb,'ylabel');
set(coll,'String','vp in m/s'); 
caxis([min_vp max_vp])
axis equal tight

subplot(3,2,4)
imagesc(X,Y,true_rho)
hold on
line(x_line,y_line,'LineStyle','--','Color','k')
plot(sources,0.2*ones(1,5),'*r','MarkerSize',12)
hold off
xlabel('x in m','FontSize',fontsize)
ylabel('y in m','FontSize',fontsize)
title(['true density model in kg/m**3'],'FontSize',fontsize)
colb=colorbar;
coll=get(colb,'ylabel');
set(coll,'String','\rho in kg/m^3'); 
caxis([min_rho max_rho])
axis equal tight

% Plot initial model
subplot(3,2,3)
imagesc(X,Y,start_vs)
hold on
line(x_line,y_line,'LineStyle','--','Color','k')
plot(sources,0.2*ones(1,5),'*r','MarkerSize',12)
hold off
xlabel('x in m','FontSize',fontsize)
ylabel('y in m','FontSize',fontsize)
title(['initial S-wave velocity model in m/s'],'FontSize',fontsize)
colb=colorbar;
coll=get(colb,'ylabel');
set(coll,'String','vs in m/s'); 
caxis([min_vs max_vs])
axis equal tight

% imagesc(X,Y,start_vp)
% hold on
% line(x_line,y_line,'LineStyle','--','Color','k')
% plot(sources,0.2*ones(1,5),'*r','MarkerSize',12)
% hold off
% xlabel('x in m','FontSize',fontsize)
% ylabel('y in m','FontSize',fontsize)
% title(['initial P-wave velocity model in m/s'],'FontSize',fontsize)
% colb=colorbar;
% coll=get(colb,'ylabel');
% set(coll,'String','vp in m/s'); 
% caxis([min_vp max_vp])
% axis equal tight

% imagesc(X,Y,start_rho)
% hold on
% line(x_line,y_line,'LineStyle','--','Color','k')
% plot(sources,0.2*ones(1,5),'*r','MarkerSize',12)
% hold off
% xlabel('x in m','FontSize',fontsize)
% ylabel('y in m','FontSize',fontsize)
% title(['initial density model in kg/m**3'],'FontSize',fontsize)
% colb=colorbar;
% coll=get(colb,'ylabel');
% set(coll,'String','\rho in kg/m^3'); 
% caxis([min_rho max_rho])
% axis equal tight

% Plot obtained model
subplot(3,2,5)
imagesc(X,Y,inv_vs)
hold on
line(x_line,y_line,'LineStyle','--','Color','k')
plot(sources,0.2*ones(1,5),'*r','MarkerSize',12)
hold off
xlabel('x in m','FontSize',fontsize)
ylabel('y in m','FontSize',fontsize)
title(['obtained S-wave velocity model in m/s (iteration step ' int2str(iteration_no) ')'],'FontSize',fontsize)
colb=colorbar;
coll=get(colb,'ylabel');
set(coll,'String','vs in m/s'); 
caxis([min_vs max_vs])
axis equal tight

% Plot vertical velocity provfiles for S-wave velocity model
x1=20;
x2=25;

subplot(3,2,6)
plot(true_vs(:,round(x1/DH)),Y,'Color',[0.7 0.7 0.7],'LineWidth',6,'LineStyle','-');
hold on
plot(start_vs(:,round(x1/DH)),Y,'k','LineWidth',3,'LineStyle','--');
plot(inv_vs(:,round(x1/DH)),Y,'b','LineWidth',3);
plot(inv_vs(:,round(x2/DH)),Y,'r','LineWidth',3,'LineStyle','--');
line([min_vs-10 max_vs+10],[(ny-FW)*DH (ny-FW)*DH],'LineStyle','--','Color','k')
hold off
xlim([min_vs-10 max_vs+10])
set(gca,'ydir','reverse')
xlabel('velocity in m/s','fontsize',fontsize)
ylabel('depth in m','fontsize',fontsize)
title('Vertical velocity profiles','fontsize',fontsize)
legend('true','initial',['x=' num2str(x1) 'm'],['x=' num2str(x2) 'm']);
hold off
grid on
