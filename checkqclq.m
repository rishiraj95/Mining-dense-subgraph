function check = checkqclq(Y,G,gamma)

GY=subgraph(G,Y);
    for i=1:length(Y)
        if degree(GY,i)<gamma*(length(Y)-1)
            check=false;
            return
        end
     
    end
    
    check=true;
    
end