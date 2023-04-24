function [y, label, U, iter, obj] = FCAG(B, U0)
% Input
% B   anchor graph n*m
% U0  is the initial label matrix m*c

% Output
% y     is the label vector of n samples n*1 
% label is the label vector of m anchors m*1 
% U     is the label matrix of m anchors m*c 
% iter  is the number of iteration
% obj   is the objective function value
% The code is written by Jingjing Xue 


[n,m] = size(B);
U = U0;
%% store once
aa=sum(U,1);
[~,label]=max(U,[],2);
BBB=2*(B'*B);
XX=diag(BBB)./2;

BBUU= BBB* U;% BBUU(i,:)   bh'Bul
ybby=diag(U'*BBUU/2);
%% compute Initial objective function value
  [~,Label0]=max(B*U0,[],2);
  obj(1)= sum((ybby).^(1/2));
%%
for iter=1:30  
     for i = 1:m   
         mm = label(i) ;
        if aa(mm)==1
            continue;  
        end    
       %% 
        V2=ybby'+(BBUU(i,:)+XX(i)).*(1-U(i,:));
        V1=(ybby'-(BBUU(i,:)-XX(i)).*U(i,:));
        delta= V2.^(1/2)-V1.^(1/2);  
        [~,q] = max(delta);     
        if mm~=q        
             aa(q)= aa(q) +1; %  YY(p,p)=Y(:,p)'*Y(:,p);
             aa(mm)= aa(mm) -1; %  YY(m,m)=Y(:,m)'*Y(:,m)
             ybby(mm)=V1(mm); % (30)
             ybby(q)=V2(q);
             U(i,mm)=0;
             U(i,q)=1;
             label(i)=q;
             BBUU(:,mm)=BBUU(:,mm)-BBB(:,i);% (29)
             BBUU(:,q)=BBUU(:,q)+BBB(:,i);     
        end  
     end   
%% compute objective function value
    obj(iter+1) = sum(ybby.^(1/2)) ;  
    if (iter>1 && abs(obj(iter)-obj(iter-1)) < 10^-6)
        break;
    end  
end 
    F=B*U;
    [~,y]=max(F,[],2);
end

