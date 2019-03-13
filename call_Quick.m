A=ones(10,10);
G=graph(A);
minsize=floor(size(A,1)/5);
gamma=0.4;

X=[];
deg=degree(G);

candX=find(deg>=(gamma*(minsize-1)));

[result,check]=Quick(G,X,candX,gamma,minsize);