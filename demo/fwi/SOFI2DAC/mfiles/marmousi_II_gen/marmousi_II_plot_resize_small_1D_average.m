close all
clear all

% This m-file is for making movies of wave propagation.

%%%%%%%%%%%%%%%%%%%%%%%% INPUT PARAMETER %%%%%%%%%%%%%%%%%%%%%%%%%
% Here the basic name of the binary snapshot files must be given:
% (The default extension is *.bin. The increasing number of the snapshot
% is added to the basic name, e.g. if ten snapshots were computed
% and the basic name is snap, then the filenames for the snapshots
% are snap1.bin, snap2.bin, ..., snap10.bin.
%filediv='snap/Khang_hom.bin.p.00';

% Receiver Positions
%yrec=3.6:0.66:69.6;
%xrec=45.5.*ones(length(yrec),1);
clear all
close all

% plot distance of sources and receivers
dnrec=1;
dnsour=1;

zeropad=23;     % size of the padding layer for the Gaussian filter 
water_depth=29; % water depth in grid points
filter_depth=40; % depth in which the Gaussian filter is applied 
maxfilt=100;

% write model output files
WRITEMODE=1;

rec=load('receiver_resize.dat');
xrec=rec(:,1);
yrec=rec(:,2);

xrec = xrec(1:dnrec:length(xrec));
yrec = yrec(1:dnrec:length(yrec));

source=load('sources_plot_resize.dat');
xshot=source(:,1);
yshot=source(:,3);

xshot = xshot(1:dnsour:length(xshot));
yshot = yshot(1:dnsour:length(yshot));

% -------------------------------------------------------------------------
% P-Wave Velocity
% -------------------------------------------------------------------------

figure;

% gridsize and grid spacing (as specified in parameter-file) 
NX1=1; NX2=13601;
NY1=1; NY2=2801; 
IDX=1; IDY=1;
dh=1.25;

SMOOTH=1;
SHOW=1;
IMP=1;

% time increment for snapshots:
TSNAPINC=2.0e-2; TSNAP1=0.001;
FW=0.0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% grid size
nx=NX2-NX1+1; ny=NY2-NY1+1;

% plot range and increment
xp1=NX1*dh; xp2=NX2*dh; yp1=NY1*dh; yp2=NY2*dh; 

% Computing range for axis and subscript range for the movie
x=xp1:IDX*dh:xp2;
y=yp1:IDY*dh:yp2;

% clip values for Pressure P 
%caxis_value_p_1=2178.0;
%caxis_value_p=1100.0;

% clip values for vx
%caxis_value_p_1=3000.0;
%caxis_value_p=1000.0;

%caxis_value_p_1=2350.0;
%caxis_value_p=900.0;

       
% source position
%xs=2610.0;
%ys=2100.0;

% receiver position
%xrec(1)=2610.0;
%xrec(2)=3730.0;
%yrec(1)=3230.0;
%yrec(2)=2100.0;

%firstframe=50;
%lastframe=150;

%load 'seismic.map'
colormap(jet(256));

% load model
 file='Model_vp.bin';
 disp([' loading file ' file]);
 fid=fopen(file,'r','ieee-le');
 vp=fread(fid,[ny,nx],'float');
 fclose(fid);
 
 % cut water column
 vp_wo_water = vp(246:2800,:);
 
 clear vp;
 
 size(vp_wo_water)
 
 % resample model
 vp_resamp = vp_wo_water(1:2:2475,1:2:13600);
 
 clear vp_wo_water;
 clear x;
 clear y;

vp_resamp_1 = vp_resamp(1:2:1160,2801:2:4800);
clear vp_resamp;
vp_resamp = vp_resamp_1;

size(vp_resamp)

caxis_value_vp1=1500.0;
caxis_value_vp2=4700.0;

caxis_value_vs1=0.0;
caxis_value_vs2=2713.5;

caxis_value_rho1=1000.0;
caxis_value_rho2=2566.8;
 
% gridsize and grid spacing (as specified in parameter-file) 
NX1=1; NX2=1000;
NY1=1; NY2=580; 
IDX=1; IDY=1;
dh=5.0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% grid size
nx=NX2-NX1+1 
ny=NY2-NY1+1

% plot range and increment
xp1=NX1*dh; xp2=NX2*dh; yp1=NY1*dh; yp2=NY2*dh; 

% Computing range for axis and subscript range for the movie
x=xp1:IDX*dh:xp2;
y=yp1:IDY*dh:yp2;


if(SMOOTH==1)
    
    for j=1:NY2
        avg(j) = 0.0;
        for i=1:NX2
            
            avg(j) = avg(j) + vp_resamp(j,i);
        
        end
    end

    avg = avg./NX2;    
    
% function fitting
   
   A = [y(filter_depth:length(y))' ones(1,length(y)-filter_depth+1)'];
   size(A)
   
   b=inv(A'*A)*A'*avg(filter_depth:length(y))' 
   
   avg_lin = b(1).*y(filter_depth:length(y))+b(2);
   
   avg_lin1 = [avg(1:filter_depth-1) avg_lin];
   
   for j=1:NY2
        for i=1:NX2
            
            vp_resamp(j,i)=avg_lin1(j);
        
        end
   end
    
clear avg_lin;
end



	imagesc(x,y,vp_resamp);
	hold on;
	%contour(x,y,model1,10,'k-');
    %plot(xrec,yrec,'wo');
	%plot(xshot,yshot,'r*');
	
    caxis([caxis_value_vp1 caxis_value_vp2]);
    
    colorbar;
	%set(gca,'YDir','normal');
        %colorbar;
        %axis equal;
       %set(gca,'DataAspectRatio',[1 1 1]);
       set(get(gca,'title'),'FontSize',12);
       set(get(gca,'title'),'FontWeight','bold');
       set(get(gca,'Ylabel'),'FontSize',12);
       set(get(gca,'Ylabel'),'FontWeight','bold');
       set(get(gca,'Xlabel'),'FontSize',12);
       set(get(gca,'Xlabel'),'FontWeight','bold');
       set(gca,'FontSize',12);
       set(gca,'FontWeight','bold');
       set(gca,'Box','on');
       set(gca,'Linewidth',1.0);
       axis([min(x)+FW max(x)-FW min(y) max(y)-FW])
		 axis ij
       
       xlabel('x [m]');
       ylabel('y [m]');
       
       iter_text=['P-Wave Velocity [m/s]'];
       
       title(iter_text);
   
   % function fitting
    figure;
    plot(y,avg,'b-',y,avg_lin1,'r-');
	hold on;
    
    legend('average value','linear approx.',2);
    
    set(get(gca,'title'),'FontSize',12);
       set(get(gca,'title'),'FontWeight','bold');
       set(get(gca,'Ylabel'),'FontSize',12);
       set(get(gca,'Ylabel'),'FontWeight','bold');
       set(get(gca,'Xlabel'),'FontSize',12);
       set(get(gca,'Xlabel'),'FontWeight','bold');
       set(gca,'FontSize',12);
       set(gca,'FontWeight','bold');
       set(gca,'Box','on');
       set(gca,'Linewidth',1.0);
    	 axis ij
       
       ylabel('P-wave velocity [m/s]');
       xlabel('depth [m]');
       
       iter_text=['P-Wave Velocity [m/s]'];
       
       title(iter_text);

if(WRITEMODE==1)       
file1=['marmousi_II_scaled_lin_1D.vp'];
fid1=fopen(file1,'w','ieee-le');
fwrite(fid1,vp_resamp,'float')
fclose(fid1);
end

clear x;
clear y;

% -------------------------------------------------------------------------
% S-Wave Velocity
% -------------------------------------------------------------------------
 
figure;
% gridsize and grid spacing (as specified in parameter-file) 
NX1=1; NX2=13601;
NY1=1; NY2=2801; 
IDX=1; IDY=1;
dh=1.25;


% time increment for snapshots:
TSNAPINC=2.0e-2; TSNAP1=0.001;
FW=0.0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% grid size
nx=NX2-NX1+1; ny=NY2-NY1+1;

% plot range and increment
xp1=NX1*dh; xp2=NX2*dh; yp1=NY1*dh; yp2=NY2*dh; 

% Computing range for axis and subscript range for the movie
x=xp1:IDX*dh:xp2;
y=yp1:IDY*dh:yp2;

%load 'seismic.map'
colormap(jet(256));

% load model
 file='Model_vs.bin';
 disp([' loading file ' file]);
 fid=fopen(file,'r','ieee-le');
 vs=fread(fid,[ny,nx],'float');
 fclose(fid);
 
 % cut water column
 vs_wo_water = vs(246:2800,:);
 vs_wo_water(1:160,:) = 303.4;
 clear vs;
 
 size(vs_wo_water)
 
 % resample model
 vs_resamp = vs_wo_water(1:2:2475,1:2:13600); 
 
 clear vs_wo_water;
 clear x;
 clear y;

vs_resamp_1 = vs_resamp(1:2:1160,2801:2:4800);
clear vs_resamp;
vs_resamp = vs_resamp_1;
 
 size(vs_resamp) 
 
% gridsize and grid spacing (as specified in parameter-file) 
NX1=1; NX2=1000;
NY1=1; NY2=580; 
IDX=1; IDY=1;
dh=5.0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% grid size
nx=NX2-NX1+1; ny=NY2-NY1+1;

% plot range and increment
xp1=NX1*dh; xp2=NX2*dh; yp1=NY1*dh; yp2=NY2*dh; 

% Computing range for axis and subscript range for the movie
x=xp1:IDX*dh:xp2;
y=yp1:IDY*dh:yp2;

vs_resamp = vp_resamp ./ sqrt(3);

vs_resamp(1:water_depth,:) = 0.0;

if(SMOOTH==1)
    
    for j=1:NY2
        avg(j) = 0.0;
        for i=1:NX2
            
            avg(j) = avg(j) + vs_resamp(j,i);
        
        end
    end
    
    avg = avg./NX2;    

   figure;
   % function fitting
   
   A = [y(filter_depth:length(y))' ones(1,length(y)-filter_depth+1)'];
   size(A)
   
   b=inv(A'*A)*A'*avg(filter_depth:length(y))' 
   
   avg_lin = b(1).*y(filter_depth:length(y))+b(2);
   
   avg_lin1 = [avg(1:filter_depth-1) avg_lin];
   
   for j=1:NY2
        for i=1:NX2
            
            vs_resamp(j,i)=avg_lin1(j);
        
        end
    end
    
end



	imagesc(x,y,vs_resamp);
	hold on;
	%contour(x,y,model1,10,'k-');
    %plot(xrec,yrec,'wo');
	%plot(xshot,yshot,'w*');
	
    caxis([caxis_value_vs1 caxis_value_vs2]);
    
    colorbar;
	%set(gca,'YDir','normal');
        %colorbar;
        %axis equal;
       %set(gca,'DataAspectRatio',[1 1 1]);
       set(get(gca,'title'),'FontSize',12);
       set(get(gca,'title'),'FontWeight','bold');
       set(get(gca,'Ylabel'),'FontSize',12);
       set(get(gca,'Ylabel'),'FontWeight','bold');
       set(get(gca,'Xlabel'),'FontSize',12);
       set(get(gca,'Xlabel'),'FontWeight','bold');
       set(gca,'FontSize',12);
       set(gca,'FontWeight','bold');
       set(gca,'Box','on');
       set(gca,'Linewidth',1.0);
       axis([min(x)+FW max(x)-FW min(y) max(y)-FW])
		 axis ij
       
       xlabel('x [m]');
       ylabel('y [m]');
       
       iter_text=['S-Wave Velocity [m/s]'];
       
       title(iter_text);
       
       hold off;

    figure;
    plot(y,avg,'b-',y,avg_lin1,'r-');
	hold on;
    
    legend('average value','linear approx.',2);
    hold on;

    set(get(gca,'title'),'FontSize',12);
       set(get(gca,'title'),'FontWeight','bold');
       set(get(gca,'Ylabel'),'FontSize',12);
       set(get(gca,'Ylabel'),'FontWeight','bold');
       set(get(gca,'Xlabel'),'FontSize',12);
       set(get(gca,'Xlabel'),'FontWeight','bold');
       set(gca,'FontSize',12);
       set(gca,'FontWeight','bold');
       set(gca,'Box','on');
       set(gca,'Linewidth',1.0);
    	 axis ij
       
       ylabel('S-wave velocity [m/s]');
       xlabel('depth [m]');
       
       iter_text=['S-Wave Velocity [m/s]'];
       
       title(iter_text);

if(WRITEMODE==1) 
file1=['marmousi_II_scaled_lin_1D.vs'];
fid1=fopen(file1,'w','ieee-le');
fwrite(fid1,vs_resamp,'float')
fclose(fid1);
end

clear x;
clear y;
clear vs_resamp;

% -------------------------------------------------------------------------
% Density
% -------------------------------------------------------------------------

% gridsize and grid spacing (as specified in parameter-file) 
NX1=1; NX2=13601;
NY1=1; NY2=2801; 
IDX=1; IDY=1;
dh=1.25;

FW=0.0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% grid size
nx=NX2-NX1+1; ny=NY2-NY1+1;

% plot range and increment
xp1=NX1*dh; xp2=NX2*dh; yp1=NY1*dh; yp2=NY2*dh; 

% Computing range for axis and subscript range for the movie
x=xp1:IDX*dh:xp2;
y=yp1:IDY*dh:yp2;

%load 'seismic.map'
colormap(jet(256));

% load model
 file='Model_rho.bin';
 disp([' loading file ' file]);
 fid=fopen(file,'r','ieee-le');
 rho=fread(fid,[ny,nx],'float');
 fclose(fid);
 
 % cut water column
 rho_wo_water = rho(246:2800,:);
 rho_wo_water(1:160,:) = 1.956;
 
 clear rho;
 
 size(rho_wo_water)
 
 % resample model
 rho_resamp = rho_wo_water(1:2:2475,1:2:13600).*1000.0; 
 
 clear rho_wo_water;
 clear x;
 clear y;
 
rho_resamp_1 = rho_resamp(1:2:1160,2801:2:4800);
clear rho_resamp;
rho_resamp = rho_resamp_1;
 
 size(rho_resamp)

 
% gridsize and grid spacing (as specified in parameter-file) 
NX1=1; NX2=1000;
NY1=1; NY2=580; 
IDX=1; IDY=1;
dh=5.0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% grid size
nx=NX2-NX1+1; ny=NY2-NY1+1;

% plot range and increment
xp1=NX1*dh; xp2=NX2*dh; yp1=NY1*dh; yp2=NY2*dh; 

% Computing range for axis and subscript range for the movie
x=xp1:IDX*dh:xp2;
y=yp1:IDY*dh:yp2;

rho_resamp = 1000.0 .* 0.31 .* (vp_resamp).^(1./4);

rho_resamp(1:water_depth,:) = 1000.0;


if(SMOOTH==1)
    
    for j=1:NY2
        avg(j) = 0.0;
        for i=1:NX2
            
            avg(j) = avg(j) + rho_resamp(j,i);
        
        end
    end
    
    avg = avg./NX2;    
    
    % function fitting
   
   A = [y(filter_depth:length(y))' ones(1,length(y)-filter_depth+1)'];
   size(A)
   
   b=inv(A'*A)*A'*avg(filter_depth:length(y))' 
   
   avg_lin = b(1).*y(filter_depth:length(y))+b(2);
   
   avg_lin1 = [avg(1:filter_depth-1) avg_lin];
   
   for j=1:NY2
        for i=1:NX2
            
            rho_resamp(j,i)=avg_lin1(j);
        
        end
   end

end

figure; 
	imagesc(x,y,rho_resamp);
    
    hold on;
	%contour(x,y,model1,10,'k-');
    %plot(xrec,yrec,'wo');
	%plot(xshot,yshot,'w*');
	
    caxis([caxis_value_rho1 caxis_value_rho2]);
    
    colorbar;
	%set(gca,'YDir','normal');
        %colorbar;
        %axis equal;
       %set(gca,'DataAspectRatio',[1 1 1]);
       set(get(gca,'title'),'FontSize',12);
       set(get(gca,'title'),'FontWeight','bold');
       set(get(gca,'Ylabel'),'FontSize',12);
       set(get(gca,'Ylabel'),'FontWeight','bold');
       set(get(gca,'Xlabel'),'FontSize',12);
       set(get(gca,'Xlabel'),'FontWeight','bold');
       set(gca,'FontSize',12);
       set(gca,'FontWeight','bold');
       set(gca,'Box','on');
       set(gca,'Linewidth',1.0);
       axis([min(x)+FW max(x)-FW min(y) max(y)-FW])
		 axis ij
       
       xlabel('x [m]');
       ylabel('y [m]');
       
       iter_text=['Density [kg/m^3]'];
       
       title(iter_text);

    figure;
    plot(y,avg,'b-',y,avg_lin1,'r-');
	hold on;
    
    legend('average value','linear approx.',2);
    
   
	hold on;

    set(get(gca,'title'),'FontSize',12);
       set(get(gca,'title'),'FontWeight','bold');
       set(get(gca,'Ylabel'),'FontSize',12);
       set(get(gca,'Ylabel'),'FontWeight','bold');
       set(get(gca,'Xlabel'),'FontSize',12);
       set(get(gca,'Xlabel'),'FontWeight','bold');
       set(gca,'FontSize',12);
       set(gca,'FontWeight','bold');
       set(gca,'Box','on');
       set(gca,'Linewidth',1.0);
    	 axis ij
       
       ylabel('Density [kg/m^3]');
       xlabel('depth [m]');
       
       iter_text=['Density [kg/m^3]'];
       
       title(iter_text);
       
       
if(WRITEMODE==1)       
% rho_resamp(1,:) = 1.25;
% rho_resamp(2,:) = 1.25;

file1=['marmousi_II_scaled_lin_1D.rho'];
fid1=fopen(file1,'w','ieee-le');
fwrite(fid1,rho_resamp,'float')
fclose(fid1);
end




















