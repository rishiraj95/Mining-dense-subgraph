function [result, ismaxqclq] = Quick(G,X,candX,gamma,minsize)

    ismaxqclq=false;
    result={};
    knear=4;
    Cu_final=[];
    %Find cover vertex set
    for i=1:length(candX)
        u=candX(i);
        notNu=setdiff((1:numnodes(G)),neighbors(G,u));
        notNuX=[];
        for j=1:length(notNu)
            if (ismember(notNu(j),X))
                notNuX=[notNuX,notNu(j)];
            end
        end

        NnotNuX=[];
        for j=1:length(notNuX)
            NnotNuX=intersect(NnotNuX,neighbors(G,notNuX(j)));
        end

        Cu_temp=intersect(intersect(neighbors(G,u),candX),NnotNuX);
        if length(Cu_temp)>length(Cu_final)
            Cu_final=Cu_temp;
        end
    end

    %Put Cu_final elements at the ened of candX
    candX_nocover=setdiff(candX,Cu_final);
    candX=[candX_nocover,Cu_final];

    for i=1:length(candX_nocover)
        new_result={};
        if length(X)+length(candX)<minsize
            ismaxqclq=false;
            return
        end
        v=candX_nocover(i);
        X_union_candX=union(X,candX);
        if checkqclq(X_union_candX,G,gamma)
            new_result={X_union_candX};
%             disp(['v' num2str(v)])
%                 disp('result')
%                 for k=1:length(new_result)
%                     disp(new_result{k})
%                end
            result(length(result)+1)=new_result;
            ismaxqclq=true;
            return
        end

        Y=union(X,v);
        candX=setdiff(candX,v);
        candY=intersect(candX,nearest(G,v,knear));

        if ~isempty(candY)
            myupdate
        else
            LY=0;UY=0;
        end
        Z=1;
        while LY<=UY && ~isempty(candY) && ~isempty(Z)
            
            for j=1:length(Y)
                indeg=indeg_Y(j);
                exdeg=exdeg_Y(j);

                if (indeg+exdeg)==gamma*(length(Y)+LY-1)
                    addnow=intersect(candY,neighbors(G,Y(j)));
                    Y=union(Y,addnow);
                    candY=setdiff(candY,addnow);
                    if ~isempty(candY)
                         myupdate
                    else
                        LY=0;UY=0;
                    end
                    break
                end
            end

          

            Z=Y((indeg_Y+exdeg_Y<gamma*(length(Y)+exdeg_Y-1))|(indeg_Y+UY<gamma*(length(Y)+UY-1))|(indeg_Y+exdeg_Y<gamma*(length(Y)+LY-1)));
            if ~isempty(Z)
                candY=[];
                
            else
                Z=candY((indeg_candY+exdeg_candY<gamma*(length(Y)+exdeg_candY))|(indeg_candY+UY-1<gamma*(length(Y)+UY-1))|(indeg_candY+exdeg_candY<gamma*(length(Y)+LY-1)));
                candY=setdiff(candY,Z);
            end
            
            if ~isempty(candY)
                myupdate
            else
                LY=0;UY=0;
            end
        end


        if LY<=UY && ~isempty(candY) && (length(Y)+length(candY))>minsize
            [new_result, issuperqclq]=Quick(G,Y,candY,gamma,minsize);
%             disp(['v' num2str(v)])
%                 disp('new_result')
%                 for k=1:length(new_result)
%                     disp(new_result{k})
%                 end
            ismaxqclq=ismaxqclq | issuperqclq;
            if length(Y)>minsize && checkqclq(Y,G,gamma) && ~issuperqclq
                ismaxqclq=true;
                new_result={Y};
%                 disp(['v' num2str(v)])
%                 disp('new_result')
%                 for k=1:length(new_result)
%                     disp(new_result{k})
%                 end
            end
            
        end
        
%         disp(['v' num2str(v)])
%         disp('result')
%                 for k=1:length(result)
%                     disp(result{k})
%                 end 
        
        
            
            if ~isempty(new_result)
                result(length(result)+1:length(result)+length(new_result))=new_result;
%                 disp(['v' num2str(v)])
%                 disp('result')
%                 for k=1:length(new_result)
%                     disp(new_result{k})
%                 end
                    
            end
        
    end
end
            
        
        
        
        
            
            
        
          
        



    
    



