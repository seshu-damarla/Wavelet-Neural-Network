
%dvenugopalarao%

function ypred=wnn(n,wi,a,b,wo,xnew)

m=1;
Fb=0.5;
Fc=0.5;

net=zeros(1,n);
net_ab=zeros(1,n);

M=size(xnew,2);

Q=size(xnew,1);
ypred=zeros(Q,1);

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
        
     output=wo*H';
     ypred(i,1)=output;
     
     net=zeros(1,n);
     net_ab=zeros(1,n);
 
end