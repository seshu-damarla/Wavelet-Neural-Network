%dvgpro%

function [Xtraining_data,Xtesting_data,training_QualityVar,testing_Qualityvar]=train_test_data(x,y,method,Ns,sigma)
        
% the data is divided into trainig set, and testing set
Y_actual=y;X_actual=x;
if strcmp(method,'RS')
    
    idx=randperm(length(Y_actual));
    size1=Ns;
    indexToTraing=(idx<=size1);
    indexToTesting=(idx>size1);
    Xtraining_data=X_actual(indexToTraing,:); % input to SFA
    Xtesting_data=X_actual(indexToTesting,:); % input to SFA

    training_QualityVar=Y_actual(indexToTraing,:);
    testing_Qualityvar=Y_actual(indexToTesting,:);

elseif strcmp(method,'KS')
    
    [model,test]=kenstone(X_actual,Ns);
    
    Xtraining_data=X_actual(model,:);
    Xtesting_data=X_actual(test,:);
    
    training_QualityVar=Y_actual(model,:);
    testing_Qualityvar=Y_actual(test,:);
    
elseif strcmp(method,'DS')
    [model,test] = Duplex(X_actual,Ns);
    Xtraining_data=X_actual(model,:);
    Xtesting_data=X_actual(test,:);
    
    training_QualityVar=Y_actual(model,:);
    testing_Qualityvar=Y_actual(test,:);
    
elseif strcmp(method,'ES')
    
    Xtraining_data=X_actual(1:Ns,:); % input to SFA
    Xtesting_data=X_actual(Ns+1:end,:); % input to SFA

    training_QualityVar=Y_actual(1:Ns,:);
    testing_Qualityvar=Y_actual(Ns+1:end,:);
    
elseif strcmp(method,'HS')
    
    [model,test] = HSPXY(X_actual,Y_actual,Ns);
    Xtraining_data=X_actual(model,:);
    Xtesting_data=X_actual(test,:);
    
    training_QualityVar=Y_actual(model,:);
    testing_Qualityvar=Y_actual(test,:);
    
elseif strcmp(method,'SP')
    [model,test] = spxy(X_actual,Y_actual,Ns);
    Xtraining_data=X_actual(model,:);
    Xtesting_data=X_actual(test,:);
    
    training_QualityVar=Y_actual(model,:);
    testing_Qualityvar=Y_actual(test,:);   
    
elseif strcmp(method,'KSP')
    [model,test] = kspxy(X_actual,Y_actual,Ns,sigma);
    Xtraining_data=X_actual(model,:);
    Xtesting_data=X_actual(test,:);
    
    training_QualityVar=Y_actual(model,:);
    testing_Qualityvar=Y_actual(test,:);    
end

% if Nm==1
%     training_QualityVar=(training_QualityVar-mean(training_QualityVar))/std(training_QualityVar);
%     testing_Qualityvar=(testing_Qualityvar-mean(testing_Qualityvar))/std(testing_Qualityvar);
%     
%     X1_actual=zeros(size(Xtraining_data));X2_actual=zeros(size(Xtesting_data));
%     meanx1=mean(Xtraining_data);stdx1=std(Xtraining_data);
%     for i=1:1:size(Xtraining_data,2)
%         X1_actual(:,i)=(Xtraining_data(:,i)-meanx1(i))/stdx1(i);
%     end
%     Xtraining_data=X1_actual;
%     meanx2=mean(Xtesting_data);stdx2=std(Xtesting_data);
%     for i=1:1:size(Xtesting_data,2)
%         X2_actual(:,i)=(Xtesting_data(:,i)-meanx2(i))/stdx2(i);
%     end
%     Xtesting_data=X2_actual;
% end

end