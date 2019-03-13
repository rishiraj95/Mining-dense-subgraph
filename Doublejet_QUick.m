clear all
close all

tout=50;
lbin=6;
numps=ncread(fullfile('../output_101.nc'),'num_particles');
numps=double(numps);
cum_adjmat=sparse(numps,numps);

for ii=1:tout
load(fullfile('../',['adjmat_' num2str(ii) '.mat']));
cum_adjmat=double(cum_adjmat|myadj);
end

G=graph(cum_adjmat);
bins=conncomp(G);
startbinsi=cell(lbin,1);

maxbins=max(bins);
for i=1:maxbins
    binpop(i)=sum(bins==i);
end
[binsort,binsorti]=sort(binpop,'descend');

for topi=1:lbin
 binnowi=find(bins==binsorti(topi));
 startbinsi{topi}=binnowi;
end

partnow=startbinsi{4};
Gnow=subgraph(G,partnow);
disp(length(partnow));
minsize=10;
gamma=0.3;

X=[];
deg=degree(Gnow);

candX=find(deg>=(gamma*(minsize-1)));

[result,check]=Quick(Gnow,X,candX,gamma,minsize);

ii=50;
part_x=ncread(fullfile('../',['output_' num2str(ii) '.nc']),'particle_x_position');
part_y=ncread(fullfile('../',['output_' num2str(ii) '.nc']),'particle_y_position');
figure(1)
clf
for i=1:length(result)
ploti=startbinsi{2}(result{i});
plot(part_x(ploti),part_y(ploti),'bo','MarkerSize',6), axis([1 9 1 7.5]), hold on
plot(part_x(startbinsi{2}),part_y(startbinsi{2}),'k.','MarkerSize',4)
drawnow
hold off
pause(0.4)
end

