closereq;
clear all;
clc;

dd=33;                               % total points number
dox=2500;                            % whole x domain
bottom=0;                        % bottom y domain
ztop=2500;                            % up y domain

tout=5;                             % total output number
tmax=0.075*3200;                         % time max

timelabel=0:tmax/(tout-1):tmax;      % time increasing step

ww=1.4;
mx=4.5;                             % highest altitude(m)

%---------- open hh.dat -------------------------
ftr=fopen('../13/temp_asc.dat');
[rdatatr,count]=fscanf(ftr,'%e');
%data=fscanf(ftr,'%e');
a=count/tout;

for m=tout:-1:1
    tr(:,m)=rdatatr((a*(m-1)+1):(a*m),1);
    pp=tr(:,m);
    pp=reshape(pp,dd,dd);
    tr2(:,:,m)=pp';
end

%---------- open uu.dat -------------------------
ftr=fopen('../13/uu_asc.dat','r');
[rdatatr,count]=fscanf(ftr,'%e');
a=count/tout;

for m=tout:-1:1
    uu(:,:,m)=rdatatr((a*(m-1)+1):(a*m),1);
    qq=uu(:,:,m);
    qq=reshape(qq,dd,dd);
    uu2(:,:,m)=qq';
end

%---------- open ww.dat -------------------------
ftr=fopen('../13/ww_asc.dat','r');
[rdatatr,count]=fscanf(ftr,'%e');
a=count/tout;

for m=tout:-1:1
    vv(:,:,m)=rdatatr((a*(m-1)+1):(a*m),1);
    rr=vv(:,:,m);
    rr=reshape(rr,dd,dd);
    vv2(:,:,m)=rr';
end

%------------------------------------------------

status=fclose('all');

x=0:dox/(dd-1):dox;
xx=0:2*dox/(dd-1):dox;
y=bottom:(ztop-bottom)/(dd-1):ztop;
yy=bottom:2*(ztop-bottom)/(dd-1):ztop;

%colormap([0 0 1; 0.1 0.1 1; 0.2 0.2 1; 0.3 0.3 1; 0.4 0.4 1; 0.5 0.5 1; 0.6 0.6 1; 0.7 0.7 1; 0.8 0.8 1; 0.9 0.9 1; 
%          1 1 0.9; 1 1 0.8; 1 1 0.7; 1 1 0.6; 1 1 0.5; 1 1 0.4; 1 1 0.3; 1 1 0.2; 1 1 0.1; 1 1 0]);

for i=tout:-1:1
    figure(i);
    set(i,'position',[0 40 1024 650]);           % for 1024 *  768
%    set(i,'position',[0 40 1280 906]);            % for 1208 * 1024
    axes('fontsize',18);
    [c,h]=contourf(x,y,tr2(:,:,i),10);
%    set(h,'EdgeColor',[1 1 1]);            % white
%    set(h,'EdgeColor','none');
    hold on;
%    quiver(xx,yy,uu2(1:as:dd,1:as:dd,i),vv2(1:as:dd,1:as:dd,i),ww,'k');
    quiver(x,y,uu2(:,:,i),vv2(:,:,i),ww,'k');
    caxis([-mx mx]);
    colorbar('vert');
    axis square;
    axis([0,dox,bottom,ztop]);

    xlabel('X (m)','fontsize',22);
    ylabel('Z (m)','fontsize',22);

    title(['Time =   ',num2str(timelabel(i)),' seconds'],'fontsize',22);
%    output
    eval(['saveas(figure(i),''D:\plots\','tempwind',num2str(i),'.tiff'');'])
end