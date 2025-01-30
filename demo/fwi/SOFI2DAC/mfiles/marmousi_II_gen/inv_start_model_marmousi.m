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

rec=load('receiver.dat');
xrec=rec(:,1);
yrec=rec(:,2);

source=load('sources_plot.dat')
xshot=source(:,1);
yshot=source(:,3);

% gridsize and grid spacing (as specified in parameter-file) 
NX1=1; NX2=1000;
NY1=1; NY2=580; 
IDX=1; IDY=1;
dh=5.0;

SMOOTH=0;
SHOW=3;
IMP=0;

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

firstframe=100;
lastframe=100;
snapno=0;


%load 'seismic.map'
colormap(jet(256));



caxis_value_vp1=1500.0;
caxis_value_vp2=4700.0;

caxis_value_vs1=0.0;
caxis_value_vs2=2713.5;

caxis_value_rho1=1000.0;
caxis_value_rho2=2566.8;

if(IMP==1)
caxis_value_vp1=min(min(vp.*rho));
caxis_value_vp2=max(max(vp.*rho));

caxis_value_vs1=min(min(vs.*rho));
caxis_value_vs2=max(max(vs.*rho));

caxis_value_rho1=min(min(rho));
caxis_value_rho2=max(max(rho));
end

if(IMP==1)
  vp=rho.*vp;
end

if(IMP==1)
  vs=rho.*vs;
end

h1=1;
for k=firstframe:1:lastframe,

if(IMP==1)    
% load model
 file1=['waveform_test_model_rho_it_',int2str(k),'.bin'];
 disp([' loading file ' file1]);
 fid1=fopen(file1,'r','ieee-le');
 rho1=fread(fid1,[ny,nx],'float');
 fclose(fid1);
end
    
if(SHOW==1)

% load model
 file1=['model_result\waveform_test_model_vp_it_',int2str(k),'.bin'];
 disp([' loading file ' file1]);
 fid1=fopen(file1,'r','ieee-le');
 vp1=fread(fid1,[ny,nx],'float');
 fclose(fid1);

if(IMP==1)
    vp1=rho1.*vp1;
end

	imagesc(x,y,vp1);
	hold on;
    %plot(xrec,yrec,'wo');
	%plot(xshot,yshot,'w*');
	
    caxis([caxis_value_vp1 caxis_value_vp2]);
    
    %colorbar;
	%set(gca,'YDir','normal');
        colorbar;
        %axis equal;
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
       
       iter_text=['P-Wave Velocity [m/s](Iteration 100)'];
       title(iter_text);
       
       hold off;


   	pause(0.1);
   
    if(IMP==0)
    ppmfile=['pics/Vp/vp_inv_',int2str(h1),'.tiff'];
    eval(['print -dtiff ' ppmfile]);
    end
    
    if(IMP==1)
    ppmfile=['pics/Zp/Zp_inv_',int2str(h1),'.tiff'];
    eval(['print -dtiff ' ppmfile]);
    end
    
    h1=h1+1;

end

if(SHOW==2)

% load model
 file1=['model_result\waveform_test_model_vs_it_',int2str(k),'.bin'];
 disp([' loading file ' file1]);
 fid1=fopen(file1,'r','ieee-le');
 vs1=fread(fid1,[ny,nx],'float');
 fclose(fid1);

if(IMP==1)
    vs1=rho1.*vs1;
end

	imagesc(x,y,vs1);
	hold on;
	%contour(x,y,model1,10,'k-');
    %plot(xrec,yrec,'wo');
	%plot(xshot,yshot,'w*');
	caxis([caxis_value_vs1 caxis_value_vs2]);
	%colorbar;
	%set(gca,'YDir','normal');
        colorbar;
        %axis equal;
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
       
       iter_text=['S-Wave Velocity [m/s] (Iteration 100)'];
       title(iter_text);
       hold off;

   	pause(0.1);
    
    if(IMP==0)
    %ppmfile=['pics/Vs/vs_inv_',int2str(h1),'.tiff'];
    %eval(['print -dtiff ' ppmfile]);
    end
    
    if(IMP==1)
    ppmfile=['pics/Zs/Zs_inv_',int2str(h1),'.tiff'];
    eval(['print -dtiff ' ppmfile]);
    end
    
    h1=h1+1;

end

if(SHOW==3)

% load model
 file1=['model_result\waveform_test_model_rho_it_',int2str(k),'.bin'];
 disp([' loading file ' file1]);
 fid1=fopen(file1,'r','ieee-le');
 rho1=fread(fid1,[ny,nx],'float');
 fclose(fid1);

	imagesc(x,y,rho1);
	hold on;
	%contour(x,y,model1,10,'k-');
    %plot(xrec,yrec,'wo');
	%plot(xshot,yshot,'w*');
	caxis([caxis_value_rho1 caxis_value_rho2]);
	%colorbar;
	%set(gca,'YDir','normal');
        colorbar;
        %axis equal;
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
       
       iter_text=['Density [kg/m^3](Iteration 100)'];
       title(iter_text);
       hold off;

   	pause(0.1);
    
    %ppmfile=['pics/Rho/rho_inv_',int2str(h1),'.tiff'];
    %eval(['print -dtiff ' ppmfile]);
    h1=h1+1;

end
   
if(SHOW==4)
figure;       



	imagesc(x,y,taper_coeff);
	hold on;
	%contour(x,y,model1,10,'k-');
    plot(xrec,yrec,'wo');
	plot(xshot,yshot,'w*');
	%caxis([caxis_value_p caxis_value_p_1]);
	colorbar;
	%set(gca,'YDir','normal');
        %colorbar;
        %axis equal;
       set(gca,'DataAspectRatio',[1 1 1]);
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
       
       %iter_text=['Iteration step No.',int2str(i)];
       %af=text(245,-50,iter_text);
       %set(af,'FontSize',12,'FontWeight','bold');
       title('Taper coefficient []');
       
       hold off;

   	pause(0.1);
end

	%set(gcf,'Renderer','zbuffer')
	%M(i-firstframe+1)=getframe(gcf);

	%brighten(0.5)
	

	snapno=snapno+1;
   % Saving the snapshot:
%epsfile=['Khang_grid_',int2str(i),'.eps'];
%eval(['print -depsc ' epsfile]);
%    jpgfile=['jpeg/pi_it_',int2str(i),'.jpg'];
%    eval(['print -djpeg100 ' jpgfile]);
end
 
% model output 
%imfile=['crase_smooth_model_vp.dat'];
%fid = fopen(imfile,'w');

%imfile1=['crase_smooth_model_vs.dat'];
%fid1 = fopen(imfile1,'w');

%imfile2=['crase_smooth_model_rho.dat'];
%fid2 = fopen(imfile2,'w');


%for i=1:NX2
%   for j=1:NY2

%   fprintf(fid,'%e\n',vp(j,NX2./2));
%   fprintf(fid1,'%e\n',vs(j,NX2./2));
%   fprintf(fid2,'%e\n',rho(j,NX2./2));
   
%   end
%end











