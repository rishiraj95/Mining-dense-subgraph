function D2 = DISTFUN(ZI,ZJ)

Lx=10;
Ly=10;
M2=size(ZJ,1);

ZIrep=repmat(ZI,M2,1);

Zdiff=abs(ZIrep-ZJ);

Zdiff(:,1)=min(Zdiff(:,1),(Lx-Zdiff(:,1)));
Zdiff(:,2)=min(Zdiff(:,2),(Ly-Zdiff(:,2)));

D2=sqrt(Zdiff(:,1).^2+Zdiff(:,2).^2);

end



