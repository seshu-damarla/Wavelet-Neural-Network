
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

eta1=0.01;
eta2=0.001;

N=length(xtrain);
M=size(xtrain,2);

% constants in wavelet
m=1;
Fb=0.5;
Fc=0.5;

maxitrs=100;
J=zeros(1,maxitrs);

hidden_neurons=[6:2:75]; 
MSE=zeros(length(hidden_neurons),1);
% 5-fold cross validation to find optimum number of neighbors
ns=floor(length(xtrain)/5);

for wn=1:length(hidden_neurons)
    
    a=1;b=ns;
    n=hidden_neurons(wn);
    mse=zeros(5,1);
    
    for k=1:5
        
        b=k*ns;
        testx=xtrain(a:b,:);testy=ytrain(a:b);
        
        if k==1
            trainx=[xtrain(ns+1:end,:)];trainy=[ytrain(ns+1:end)];
        else
            trainx=[xtrain(1:(k-1)*ns,:);xtrain(k*ns+1:end,:)];trainy=[ytrain(1:(k-1)*ns);ytrain(k*ns+1:end)];
        end
            
        [trainx,mux,sigmax] = zscore(trainx);
        [trainy,muy,sigmay] = zscore(trainy);
        
        xnew=(testx-mux)./sigmax; % test dataset
        
        % initialization
        la=-1;ub=1;
%         rng(1)
        wi=la + (ub-la).*rand(n,M); % input layer weights
%         rng(2)
        bb=la + (ub-la).*rand(n,1);  % first parameter of wavelet
%         rng(3)
        aa=la + (ub-la).*rand(n,1);  % second parameter of wavelet
%         rng(4)
        wo=la + (ub-la).*rand(1,n);  % output layer weights
        
        [wi,aa,bb,wo]=wnntrain(maxitrs,trainx,trainy,n,wi,aa,bb,wo);
        ypred=wnn(n,wi,aa,bb,wo,xnew);
        ypred=ypred*sigmay+muy;
        
        mse(k)=mean((testy-ypred).^2);
        a=b+1;
        
        close
    end
    
    MSE(wn,1)=mean(mse);
    clear mse
    
end

plot(hidden_neurons,MSE,'LineWidth',1.5)
xlabel('Number of hidden neurons, L')
ylabel('Mean squared error')


