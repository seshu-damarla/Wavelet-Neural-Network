
%dvenugopalarao%

function [wi,a,b,wo]=wnntrain(maxitrs,xtrainn,ytrainn,n,wi,a,b,wo)

m=1;
Fb=0.5;
Fc=0.5;

eta1=0.01;
eta2=0.001;

N=length(xtrainn);
M=size(xtrainn,2);

dwi=zeros(n,M);
dwo=zeros(1,n);
db=zeros(n,1);
da=zeros(n,1);

net=zeros(1,n);
net_ab=zeros(1,n);
H=zeros(n,1); 
H1=zeros(n,1);
J=zeros(1,maxitrs);

h=animatedline;    % for dynamic plot
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
