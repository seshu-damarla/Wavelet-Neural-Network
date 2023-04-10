
%dvenugopalarao%

function J=LSfn(x)

global xtrainn ytrainn M n N Fb Fc m error count h

% M -- number of input variables
% n -- number of hidden neurons
% N -- number of training examples

xtrain=xtrainn;
ytrain=ytrainn;

wi=zeros(n,M); % input layer weights
index=1;
for i=1:n
    for j=1:M
        wi(i,j)=x(index);
        index=index+1;
    end
end
clear i j

b=zeros(n,1);  % first parameter of wavelet
for i=1:n
    b(i,1)=x(index);
    index=index+1;
end
clear i

a=zeros(n,1);  % second parameter of wavelet
for j=1:n
    a(j,1)=x(index);
    index=index+1;
end
clear j

wo=zeros(1,n);  % output layer weights
for i=1:n
    wo(i)=x(index);
    index=index+1;
end
clear i

net=zeros(1,n);
net_ab=zeros(1,n);
H=zeros(n,1);      % output of hidden layer neurons
error=zeros(N,1);        % residuals

for i=1:N
    
    inputx=xtrain(i,:);
    outputy=ytrain(i,:);
%     disp(outputy)
    
    for j=1:n 
        for k=1:M
            net(j)=net(j)+wi(j,k)*inputx(k);
        end
        net_ab(j)=(net(j)-b(j))/a(j);
        temp=fbspwavf1(net_ab(j),m,Fb,Fc);
        temp=real(temp);
        H(j)=temp;
    end
    output=wo*H;
%     disp(output)
    error(i,1)=outputy-output;
%     disp(error(i,1))
end
% error;
J=sum((error.^2)/2);  % least squares loss function

count=count+1;

addpoints(h,count,J);
drawnow

end






