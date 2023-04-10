
%dvenugopalarao%

clc
clear
close all

% rng default

[data,txt,raw]=xlsread('N2-Data.xlsx',1);
        
x=[data(:,1:4) data(:,5)];
y=data(:,7);

Ns=floor(0.8*length(data));
[xtrain,xtest,ytrain,ytest]=train_test_data(x,y,'HS',Ns,0);


M=size(xtrain,2); % no. of input variables
N=size(xtrain,1); % no.of training examples

minx=min(xtrain);maxx=max(xtrain);
xtrainn=[(xtrain-minx)./(maxx-minx)]*(1-(-1))+(-1);

miny=min(ytrain);maxy=max(ytrain);
ytrainn=[(ytrain-miny)./(maxy-miny)]*(1-(-1))+(-1);

n=29;  % no.of hidden neurons
eta1=0.01;
eta2=0.001;

% constants in wavelet
m=1;
Fb=0.5;
Fc=0.5;

% optimizing LS loss function using PSO

% nvars=3*n+n*M;
% fun=@LSfn;
% lb=-1*rand(1,nvars);
% ub=1*rand(1,nvars);
% maxitrs=min(10,10*nvars);
% options = optimoptions('particleswarm','SwarmSize',10,'HybridFcn',@fmincon,'MaxIterations',maxitrs);
% [x,fval,exitflag,output] = particleswarm(fun,nvars,lb,ub,options);

% la=-1;ub=1;
% x0=la + (ub-la).*rand([1 nvars]);
% x=fminunc(fun,rand([1 nvars]));

% x=ga(fun,nvars);

% % % initialization
% % la=-0.6;ub=0.6;
% % % rng(1)
% % wi=la + (ub-la).*rand(n,M); % input layer weights
% % % rng(2)
% % b=la + (ub-la).*rand(n,1);  % first parameter of wavelet
% % % rng(3)
% % a=la + (ub-la).*rand(n,1);  % second parameter of wavelet
% % % rng(4)
% % wo=la + (ub-la).*rand(1,n);  % output layer weights

load('25initwbest2.mat')
wi=wii;wo=woi;a=ai;b=bi;
% save('25initwbest3.mat','woi','wii','ai','bi')

dwi=zeros(n,M);
dwo=zeros(1,n);
db=zeros(n,1);
da=zeros(n,1);

net=zeros(1,n);
net_ab=zeros(1,n);
H=zeros(n,1); 
H1=zeros(n,1);

maxitrs=100;
J=zeros(1,maxitrs);

h=animatedline;    % for dynamic plot
xlabel('Number of epochs')
ylabel('SSE')
count=1;

for itr=1:maxitrs
    
    error=zeros(1,N);
    
    for i=1:N
        
        inputx=xtrainn(i,:);
        outputy=ytrainn(i);
        
        for j=1:n 
            for k=1:M
                net(j)=net(j)+wi(j,k)*inputx(k);
                net_ab(j)=(net(j)-b(j))/a(j);
            end
            
            temp=fbspwavf1(net_ab(j),m,Fb,Fc);
%             temp=mymorlet(net_ab(j));
            if isnan(temp)
%                 pause
            end
          temp=real(temp);
            H(j)=temp;
            H1(j)=real(d_fbspwavf1(net_ab(j),m,Fb,Fc));
%             H1(j)=real(d_mymorlet(net_ab(j)));
        end
        
        output=wo*H;
        error(i)=(outputy-output);
        
        for j=1:n         % gradients with output layer weights
            dwo(1,j)=-(outputy-output)*H(j);
        end
        for j=1:n        % gradients with input layer weights
            for k=1:M
                dwi(j,k)=-(outputy-output)*wo(j)*H1(j)*inputx(k)*(1/a(j));
            end
        end
        for j=1:n        % gradients with first parameter of wavelet
            db(j)=-(outputy-output)*wo(j)*H1(j)*(-1/a(j));
        end
        net=zeros(1,n);
        for j=1:n
            for k=1:M
                net(j)=net(j)+wi(j,k)*inputx(k);
            end
            da(j)=-(outputy-output)*wo(j)*H1(j)*(net(j)-b(j))*(-1/(a(j)^2));
        end
        
        wo=wo-eta1*dwo;  % updating output layer weights
        wi=wi-eta1*dwi;  % updating input layer weights
        b=b-eta2*db;     % updating frst parameter of wavelet
        a=a-eta2*da;     % updating second parameter of wavelet
        
        dwi=zeros(n,M);
        dwo=zeros(1,n);
        db=zeros(n,1);
        da=zeros(n,1);
    
        H=zeros(n,1); 
        H1=zeros(n,1);
        
        net=zeros(1,n);
        net_ab=zeros(1,n);
        
    end
    
    J(itr)=sum(error.^2);
%     disp(J(itr))
    addpoints(h,itr,J(itr));
    drawnow
        
end

% close

% predictions on test data
xnew=[(xtest-minx)./(maxx-minx)]*(1-(-1))+(-1);
% xnew=xtrainn;
Q=size(xnew,1);
ypred=zeros(Q,1);

net=zeros(1,n);
net_ab=zeros(1,n);

for i=1:Q
    
    inputx=xnew(i,:);
    
    for j=1:n 
        for k=1:M
            net(j)=net(j)+wi(j,k)*inputx(k);
        end
        net_ab(j)=(net(j)-b(j))/a(j);
        temp=fbspwavf1(net_ab(j),m,Fb,Fc);
        temp=real(temp);
        H(j)=temp;
%         H1(j)=real(d_fbspwavf1(net_ab(j),m,Fb,Fc));
     end
        
     output=wo*H;
     ypred(i,1)=output;
     
     net=zeros(1,n);
     net_ab=zeros(1,n);
 
end

ypred=0.5*(ypred*(maxy-miny)+(maxy-miny))+miny;

% ytrain=0.5*(ytrainn*(maxy-miny)+(maxy-miny))+miny;
% ytest=ytrain;

R2=(corr(ytest,ypred)^2);
fprintf('R^2= %4.4f \n',R2)

sse=sum((ypred-ytest).^2);
sst=sum((ytest-mean(ytest)).^2);
R2=1-(sse/sst);
% disp(R2)

AARD=100*mean(abs((ypred-ytest)./ytest));
fprintf('AARD= %4.4f \n',AARD)

RMSE=sqrt(mean((ypred-ytest).^2));
fprintf('RMSE= %4.4f \n',RMSE)

% average percent relative error
E=((ytest-ypred)./ytest)*100;
Er=mean(E);
fprintf('Er= %4.4f \n',Er)

Ea=mean(abs(E));
fprintf('Ea= %4.4f \n',Ea)


figure
plot(ytest)
hold on
plot(ypred)
hold off
legend('actual','model')


% p1 =      0.6267;
% p2 =       11.86;
% fx = p1*ypred + p2;
% 
% disp('-----------------')
% 
% sse=sum((ypred-ytest).^2);
% sst=sum((ytest-mean(ytest)).^2);
% R2=1-(sse/sst);
% disp(R2)
% 
% AARD=100*mean(abs((ypred-ytest)./ytest));
% disp(AARD)
% 
% RMSE=sqrt(mean((ypred-ytest).^2));
% disp(RMSE)
