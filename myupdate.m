indeg_Y=zeros(1,length(Y));
exdeg_Y=zeros(1,length(Y));
indeg_candY=zeros(1,length(candY));
exdeg_candY=zeros(1,length(candY));

for j=1:length(Y)
    indeg_Y(j)=sum(ismember(neighbors(G,Y(j)),Y));
    exdeg_Y(j)=sum(ismember(neighbors(G,Y(j)),candY));
end

for j=1:length(candY)
    indeg_candY(j)=sum(ismember(neighbors(G,candY(j)),Y));
    exdeg_candY(j)=sum(ismember(neighbors(G,candY(j)),candY));
end

[indeg_candY,sortindeg_i]=sort(indeg_candY,'descend');
candY=candY(sortindeg_i);

sum_indeg_Y=sum(indeg_Y);

t=0;
% disp('Y=') 
% disp(Y)
% disp('candY=')
% disp(candY)
 UYmin=((min(indeg_Y+exdeg_Y))/gamma)+1-length(Y);
% disp('indeg_Y')
% disp(indeg_Y)
% disp('indeg_candY')
% disp(indeg_candY)
% disp(['UYmin' num2str(UYmin)])
sum_indeg_candY_t=sum(indeg_candY(1:t+1));
while (sum_indeg_Y+sum_indeg_candY_t>=length(Y)*(gamma*(length(Y)+t-1))) && (t<=UYmin) 
    t=t+1;
   % disp(['t=' num2str(t)])
    if (t>=UYmin || t>=length(candY))
        break
    end
    sum_indeg_candY_t=sum(indeg_candY(1:t+1));
end

UY=t;

t=0;
indeg_Ymin=min(indeg_Y);

while indeg_Ymin+t<gamma*(length(Y)+t-1)
    t=t+1;
end


sum_indeg_candY_t=sum(indeg_candY(1:t));

while ((sum_indeg_Y+sum_indeg_candY_t)<length(Y)*gamma*(length(Y)+t-1)) && t<length(candY)+1
    t=t+1;
    sum_indeg_candY_t=sum(indeg_candY(1:t));
end

LY=t;
