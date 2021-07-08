clear all

nt=5150;
nr=800;
nz=807;
nx=807;
dt=0.0002;
x=(0:807-1)*1;
z=(0:807-1)*1;

% colormap('gray');
fid=fopen('seistzz.dat','rb');
vz=fread(fid,[nt,nr],'real*4');
fclose(fid);
figure(1);
imagesc(vz);
title('同相轴');
xlabel('distance（0.5mm）');
ylabel('Time/（0.02μs）')
zdz=vz(:,100);
for i=0:0;
zdz=zdz+vz(:,100-i)+vz(:,100+i);
end
zjf=zeros(5250,1);
t=(0:nt-1)*0.002e-05;
for jb=1:nt;
    jf=0;
    for jb2=1:jb;
        jf=jf+zdz(jb2)*dt;
    end
    zjf(jb)=jf;
end
figure(2)
%subplot(2,1,1)
%plot(t,zjf/dt,'-b')
%title('位移时间过程');
%xlabel('Time（s）');
%ylabel('distance（m）')
%subplot(2,1,2)
plot(t,zdz,'-b')
title('应力时间过程');
xlabel('Time（s）');
ylabel('v（m/s）')




%   colormap('gray');
% [nt,nr] 填充相应的数字
%
fid=fopen('snap_000500_tzz.dat','rb');
vzs1=fread(fid,[nz,nx],'real*4');
fclose(fid);
fid=fopen('snap_001700_tzz.dat','rb');
vzs2=fread(fid,[nz,nx],'real*4');
fclose(fid);
fid=fopen('snap_001750_tzz.dat','rb');
vzs3=fread(fid,[nz,nx],'real*4');
fclose(fid);
fid=fopen('snap_001800_tzz.dat','rb');
vzs4=fread(fid,[nz,nx],'real*4');
fclose(fid);
fid=fopen('snap_001850_tzz.dat','rb');
vzs5=fread(fid,[nz,nx],'real*4');
fclose(fid);
fid=fopen('snap_001900_tzz.dat','rb');
vzs6=fread(fid,[nz,nx],'real*4');
fclose(fid);
figure(4)
subplot(2,3,1)
imagesc(x,z,vzs1);
title('固体波场快照（0.33s）');
xlabel('distance（m）');
ylabel('distance（m）')
subplot(2,3,2)
imagesc(x,z,vzs2);
title('固体波场快照（0.34s）');
xlabel('distance（m）');
ylabel('distance（m）')
subplot(2,3,3)
imagesc(x,z,vzs3);
title('固体波场快照（0.35s）');
xlabel('distance（m）');
ylabel('distance（m）')
subplot(2,3,4)
imagesc(x,z,vzs4);
title('固体波场快照（0.36s）');
xlabel('distance（m）');
ylabel('distance（m）')
subplot(2,3,5)
imagesc(x,z,vzs5);
title('固体波场快照（0.37s）');
xlabel('distance（m）');
ylabel('distance（m）')
subplot(2,3,6)
imagesc(x,z,vzs6);
title('固体波场快照（0.38s）');
xlabel('distance（m）');
ylabel('distance（m）')

 %axis equal 




% [nz,nx] 填充相应的数字

%%
nt=5150;
nr=200;
dt=0.0002;
fid=fopen('wlet.dat','rb');
 wavlet=fread(fid,[nt,1],'real*4');
 fclose(fid);
 t=(0:nt-1)*dt;
 f=1/((t(2)-t(1))*nt)*(0:nt-1)/1000;
 vzf=abs(fft(wavlet));
 figure(5);
 subplot(2,1,1)
 plot(t(1:135),wavlet(1:135),'-b');
 title('子波波形');
xlabel('Time/s');
ylabel('Pa')
subplot(2,1,2)
plot(f(1:nt/10),vzf(1:nt/10),'-b')
title('子波频谱');
xlabel('frequency（KHz）');
ylabel('Pa')
% [nt] 填充相应的数字 nr=1
%由于震源子波是一维的，所以用plot即可
figure(100)
imagesc(x,z,vzs1);
title('固体波场快照（0.1s）');
xlabel('distance（m）');
ylabel('distance（m）')
%%


fid=fopen('wav.dat','rb');
wavsq=fread(fid,[nz,nx],'real*4');
fclose(fid);
figure(6)
imagesc(x,z,wavsq);
title('震源尺寸（105μs）');
xlabel('distance（mm）');
ylabel('distance（mm）')
%比较量区域波场快照
fid=fopen('snap_003250_tzz.dat','rb');
vzs22=fread(fid,[nz,nx],'real*4');
fclose(fid);
figure(9)
subplot(1,3,1)
imagesc(x,z,vzs22);
title('能量指向性分析（直角坐标）');
xlabel('distance（mm）');
ylabel('distance（mm）')
nd=vzs22(5:1204,5:1604);
q1=size(nd(:,1));
q2=size(nd(1,:));
nndd=zeros(1444,101);
r=(0:1444)*0.5;%单位毫米
thta=(0:0.01:1)*pi;%角度
for s=1:1200;
 for   z=1:1600;
   m=fix(ceil((s^2+(z-800)^2)^(1/2)));
   n=fix(acos((z-800)/(((z-800)^2+s^2)^(1/2)))/(0.01*pi)+1);
   nndd(m,n)=abs(nd(s,z));
end
end
q=(0:1444-1)*0.5;
p=(0:101-1)*0.01;
subplot(1,3,2)
imagesc(p,q,nndd);
title('能量指向性分析（极坐标）');
xlabel('弧度（π）');
ylabel('distance（mm）')
fenbu=nndd(807,:);
Fenbub=reshape(fenbu(2:101),4,length(fenbu(2:101))/4);
Fenbub=max(Fenbub);
p1=linspace(0,max(p(2:101)),length(fenbu(2:101))/4);
NNyzh1=length(p);
figure(222)
plot(p1+0.025,Fenbub,'-r')
axis([0 1 -inf inf]) 
title('能量分布与角度关系');
xlabel('弧度（π）');
ylabel('Amplitude(Pa)')