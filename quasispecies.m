T = 1e3;
N = 1e3;
n = 1:N;

% %X = [];
% %E = zeros(N,T);
% 
% for i=551:N
%     i
% 	a = rand(1,i);%exprnd(1,1,i);
% 	A = zeros(i,i);
% 	for k=1:i
%     	A(k,k) = a(1,k)./sum(a);
% 	end
%     
%     for j=1:T
% %        j
%         Q = rand(i,i);%exprnd(1,i,i);
% %         A = zeros(i,i);
% %         a = exprnd(1,1,i);
% 
%         Q = Q./sum(Q,2);
% %         for k=1:i
% %             Q(k,:) = Q(k,:)./sum(Q(k,:));
% % %             A(k,k) = a(1,k)./sum(a);
% %         end
% 
%         W = Q*A;
%         E1(i,j) = max(eig(W));
% 
% %        [V,D] = eig(W);
% %        [maxD,indD] = max(diag(D));
% %        X = [X,maxD];
%     end
%     
% end





% X = [];
%E = zeros(N,T);
% 
% for i=700:700
%     i
%     for j=1:T
%         j
%         Q = exprnd(1,i,i);
%         A = zeros(i,i);
%         a = exprnd(1,1,i);
% 
%         for k=1:i
%             Q(k,:) = Q(k,:)./sum(Q(k,:));
%             A(k,k) = a(1,k)./sum(a);
%         end
% 
%         W = Q*A;
% %        E(i,j) = max(eig(W));
% 
%        [V,D] = eig(W);
%        [maxD,indD] = max(diag(D));
%        X = [X,V(1:end,indD)];
%     end
%     
% end
% 



P = zeros(1,N);
Eps = 500;
Q = zeros(1,Eps);
eps = 1e-3;%exp((log(1e-7)/Eps):(log(1e-7)/Eps):log(1e-7));
nn=600;
for i=1:Eps
    m = mean(E(i,:));
%     m = mean(E(nn,:));
    count = 0;
    for j=1:T
        if(abs(E(i,j)-m)>eps)
%         if(abs(E(nn,j)-m)>eps(1,i))
            count = count+1;
        end
    end
    P(1,i) = count/T;
%     Q(1,i) = count/T;
end

hold on, semilogy(n,P,'o');
%figure, plot(n,sqrt(-log(P)),'o');
%hold on, semilogy(eps,Q,'o');




% 
% %P = zeros(1,N);
% Eps = 0:1e-6:1e-3;
% Q = zeros(1,length(Eps));
% %eps = 1e-6;%exp((log(1e-7)/Eps):(log(1e-7)/Eps):log(1e-7));
% nn=600;
% for i=1:length(Eps)
% %    m = mean(E(i,:));
%      m = mean(E(nn,:));
%      
%      count = length(find(abs(E(nn,1:end)-m)>=Eps(1,i)));
% 
%      
% %     count = 0;
% %     for j=1:T
% %         if(abs(E(i,j)-m)>eps)
% % %         if(abs(E(nn,j)-m)>eps(1,i))
% %             count = count+1;
% %         end
% %     end
% %    P(1,i) = count/T;
%      Q(1,i) = count/T;
% end
% 
% hold on, semilogy(Eps,Q,'o');
% % figure, plot(n,P,'o');
% % figure, plot(n,sqrt(-log(P)),'o');
% %hold on, semilogy(Eps,Q,'o');
