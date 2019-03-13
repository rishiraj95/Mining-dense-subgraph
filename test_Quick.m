%% This is the main time stepper for the particle workshop code.
clear all;
close all;
% Call the routine to define the toy velocity field
velfields_workshop2
U0 = 0.3;
V0 = 0.0; %Vertical background velocity
myper = 10;
% Set the basic parameter
dt=1e-4;
numsteps=500;
numouts=100;
% For Particles
numpsx=25;

numpsy=25;
numps=(numpsx+1)*(numpsy+1);

% Grids and initial position
x1d=linspace(0.05*Lx,0.95*Lx,numpsx+1);
y1d=linspace(0.05*Ly,0.95*Ly,numpsy+1);
[x0s,y0s]=meshgrid(x1d,y1d);
x0s=x0s(:);y0s=y0s(:);
x0g=reshape(x0s,numpsx+1,numpsy+1);
y0g=reshape(y0s,numpsx+1,numpsy+1);
xs1=x0s;ys1=y0s;

%Number of clusters
nbcluster=3;

t=0;
ts=zeros(numouts+1,1);
xs=zeros(numouts+1,numps);
ys=zeros(numouts+1,numps);
ts(1)=t;
xs(1,:)=x0s;ys(1,:)=y0s;
up1 = my_u1(t,xs1,ys1); %initial velocity x
vp1 = my_v1(t,xs1,ys1); %initial velocity y
tau = 100; %Drag timescale
% Initialize the degree matrix
mydistmat=sparse(numps,numps);
mydistmatcell=cell(numouts,1);
% the minimal distance
my_dist_eps=0.5*Lx/numpsx;
for ii=1:numouts
    for jj=1:numsteps
        t=t+dt; %note this means t is now effectively at step n+1
        % Symplectic Euler with exact for the second set of particles
        us1=U0+my_u1(t-dt,xs1,ys1)+my_u2(t-dt,xs1,ys1);
       % up1=up1+dt*(us1-up1)/tau;% dup/dt=(us-up)/tau
%        xs1=xs2+dt*us1; % usual step for x
        %xs1=mod(xs1+dt*us1,Lx); % usual step for x
        xs1=mod(xs1+dt*us1,Lx);
        vs1=V0+my_v1(t-dt,xs1,ys1)+my_v2(t-dt,xs1,ys1); % use updated x for calculating v
        %vp1=vp1+dt*(vs1-vp1)/tau;
        %ys1=mod(ys1+dt*vs1,Ly);
        ys1=mod(ys1+dt*vs1,Ly);
    end
    
    ts(ii)=t;
    xs(ii,:)=xs1;
    ys(ii,:)=ys1;
    
    
    %Calculate adjacency matrix at every time   
    poscoord = [xs1 ys1];
    
    idx=rangesearch(poscoord,poscoord,my_dist_eps,'Distance',@DISTFUN);%DISTFUN takes care of the periodic distance
       
    for i=1:numps
       mydistmat(i,idx{i})=1;
    end
    mydistmatcell{ii}=mydistmat;
   
    
end



lbin=4; %No. of bins we want to see

%Find connected componenents midway
startbinsi=cell(lbin,1);
cum_adjmat=zeros(numps,numps);
for ii=1:numouts
    cum_adjmat=double(cum_adjmat|full(mydistmatcell{ii}));
end

G=graph(cum_adjmat);
bins=conncomp(G);

%Marek's code to plot top bins.

maxbins=max(bins);
for i=1:maxbins
    binpop(i)=sum(bins==i);
end
[binsort,binsorti]=sort(binpop,'descend');

for topi=1:lbin
 binnowi=find(bins==binsorti(topi));
%  actualbinnows=zeroind(binnowi);
 startbinsi{topi}=binnowi;
end

figure(10)
plot(xs(50,startbinsi{1}),ys(50,startbinsi{1}),'k.','MarkerSize',3)
axis([0 10 0 10])
hold on

Gnow=subgraph(G,startbinsi{1});
minsize=floor(length(startbinsi{1})/10);
gamma=0.2;

X=[];
deg=degree(Gnow);

candX=find(deg>=(gamma*(minsize-1)));

[result,check]=Quick(Gnow,X,candX,gamma,minsize);

figure(9)
plot(xs(50,startbinsi{1}(result{1})),ys(50,startbinsi{1}(result{1})),'b.','MarkerSize',5)
hold on
plot(xs(50,startbinsi{1}(result{2})),ys(50,startbinsi{1}(result{2})),'r.','MarkerSize',5)
axis([0 10 0 10])