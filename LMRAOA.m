function [Best_FF,Best_P,Conv_curve]=LMRAOA(N,M_Iter,LB,UB,Dim,F_obj)
Best_P=zeros(1,Dim);
Best_FF=inf;
Conv_curve=zeros(1,M_Iter);
X1=initialization(N,Dim,UB,LB);
X=X1.*chaos(M_Iter);
Xnew=X;
Ffun=zeros(1,size(X,1));
Ffun_new=zeros(1,size(Xnew,1));

MOP_Max=1;
MOP_Min=0.2;
C_Iter=1;
Alpha=5;
Mu=0.499;


for i=1:size(X,1)
    Ffun(1,i)=F_obj(X(i,:));  
    AllFitness(i) = F_obj(X(i,:));
    if Ffun(1,i)<Best_FF
        Best_FF=Ffun(1,i);
        Best_P=X(i,:);
    end
end
    
    

while C_Iter<M_Iter+1  
    MOP=1-((C_Iter)^(1/Alpha)/(M_Iter)^(1/Alpha));   
    MOA=MOP_Min+C_Iter*((MOP_Max-MOP_Min)/M_Iter); 
    A = randi([1,N]);  
    B = randi([1,N]);
    C = randi([1,N]);
    oldfitness=AllFitness;
    [fitS,fitSind]=sort(AllFitness);
    Xa=X(fitSind(1),:);
    Xb=X(fitSind(2),:); 
    Xd=X(fitSind(3),:);
    Xc=X(fitSind(4),:);
    Xe=X(fitSind(5),:);
    C_pool=[Xa;Xb;Xd;Xc;Xe];
    Ceq=C_pool(randi(size(C_pool,1)),:);
    %Update the Position of solutions
    for i=1:size(X,1)   
           r1=rand();
            if (size(LB,2)==1)
                if r1>MOA
                    r2=rand();
                    if r2>0.5
                        Xnew(i,:)=Best_P(1,:)/(MOP+eps)*((UB-LB)*Mu+LB);
                    else
                        Xnew(i,:)=Best_P(1,:)*MOP*((UB-LB)*Mu+LB);
                    end
                    f1=2.02-i*(1.08/randi([1,Dim]));
                    X1(i,:)=Best_P(1,:)+rand*f1*(2*Ceq-X(A,:)-X(B,:));
                    XN1=F_obj(X1(i,:));
                    if XN1<F_obj(Xnew)
                        Xnew(i,:)=X1(i,:);
                    end
                else
                    r3=rand();
                    if r3>0.5
                        Xnew(i,:)=Best_P(1,:)-MOP*((UB-LB)*Mu+LB);
                    else
                        Xnew(i,:)=Best_P(1,:)+MOP*((UB-LB)*Mu+LB);
                    end
                    f1=2.02-i*(1.08/randi([1,Dim]));
                    X1(i,:)=X(C,:)+rand*f1*(X(A,:)-X(B,:)).*Levy(Dim);
                    XN1=F_obj(X1(i,:));
                    if XN1<F_obj(Xnew)
                        Xnew(i,:)=X1(i,:);
                    end           
            end
            
           
            if (size(LB,2)~=1)  
                r1=rand();
                if r1>MOA
                    r2=rand();
                    if r2>0.5
                        Xnew(i,:)=Best_P(1,:)/(MOP+eps)*((UB-LB)*Mu+LB);
                    else
                        Xnew(i,:)=Best_P(1,:)*MOP*((UB-LB)*Mu+LB);
                    end
                    f1=2.02-i*(1.08/randi([1,Dim]));
                    X1(i,:)=Best_P(1,:)+rand*f1*(2*Ceq-X(A,:)-X(B,:));
                    XN1=F_obj(X1(i,:));
                    if XN1<F_obj(Xnew)
                        Xnew(i,:)=X1(i,:);
                    end
                else
                    r3=rand();
                    if r3>0.5
                        Xnew(i,:)=Best_P(1,:)-MOP*((UB-LB).*chaos(M_Iter)+LB);
                    else
                        Xnew(i,:)=Best_P(1,:)+MOP*((UB-LB).*chaos(M_Iter)+LB);
                    end
                    f1=2.02-i*(1.08/randi([1,Dim]));
                    X1(i,:)=X(C,:)+rand*f1*(X(A,:)-X(B,:)).*Levy(Dim);%X(i,:)+
                    XN1=F_obj(X1(i,:));
                    if XN1<F_obj(Xnew)
                        Xnew(i,:)=X1(i,:);
                    end
                end               
            end
            
            end
        Flag_UB=Xnew(i,:)>UB; 
        Flag_LB=Xnew(i,:)<LB; 
        Xnew(i,:)=(Xnew(i,:).*(~(Flag_UB+Flag_LB)))+UB.*Flag_UB+LB.*Flag_LB;
        AllFitness(i) = F_obj(Xnew(i,:));
        
        if AllFitness(i)>=oldfitness(i) 
         k=2-2*C_Iter/M_Iter;
         D(i,:)=((min(X(i,:))+max(X(i,:)))/2)+((min(X(i,:))+max(X(i,:)))/2*k)-(X(i,:)/k);  
        if F_obj(D(i,:))<F_obj(X(i,:))
            Xnew(i,:)=D(i,:);
        end
        end
        
        Ffun_new(1,i)=F_obj(Xnew(i,:));  
        if Ffun_new(1,i)<Ffun(1,i)
            X(i,:)=Xnew(i,:);
            Ffun(1,i)=Ffun_new(1,i);
        end
        if Ffun(1,i)<Best_FF
        Best_FF=Ffun(1,i);
        Best_P=X(i,:);
        end
       AllFitness(i)=Best_FF;
    end
    
    Conv_curve(C_Iter)=Best_FF;
    
    if mod(C_Iter,50)==0
        display(['At iteration ', num2str(C_Iter), ' the best solution fitness is ', num2str(Best_FF)]);
    end
     
    C_Iter=C_Iter+1; 
   
end
end

function o=Levy(d)
beta=1.5;
sigma=(gamma(1+beta)*sin(pi*beta/2)/(gamma((1+beta)/2)*beta*2^((beta-1)/2)))^(1/beta);
u=randn(1,d)*sigma;v=randn(1,d);step=u./abs(v).^(1/beta);
o=step;
end
