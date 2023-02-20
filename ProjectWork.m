clc
clear all
close all
load("dati.mat")
dati = xlsread('NDX_price.xlsx','Sheet 1','A5:C161');
number_stock = 6;
tab_return= horzcat(ANF.TotalReturn,CATY.TotalReturn,FLR.TotalReturn,LQDT.TotalReturn,NUVA.TotalReturn,PII.TotalReturn);
Returns= tab_return/100;
% figure
% plot(Returns)
mu_riga = mean(Returns);
mu= transpose(mu_riga);
Sigma= cov(Returns);
%% 2
%% 2.1 computing markoviz standard
A = (ones(number_stock,1)')*(Sigma\mu);
B=mu'*(Sigma\mu);
C=ones(number_stock,1)'*(Sigma\ones(number_stock,1));
D=B*C-A^2;
g=(B*(Sigma\ones(number_stock,1))-A*(Sigma\mu))/D;
h=(C*(Sigma\mu)-A*(Sigma\ones(number_stock,1)))/D;
m=linspace(A/C,0.1,100);
Var_w0=@(m)C/D*(m-A/C).^2+1/C;
figure
plot(sqrt(Var_w0(m)),m)
hold on
%% 2.2 mu_exp
lambda= 0.005;
t=1:155;
mu_exp= sum(Returns.*(exp(-lambda*(155-t)))')'/sum(exp(-lambda*(155-t)));
%% 2.3 stimo matrice di cov con COSTANT CORRELATION APPROACH(CC)
Correlations=corr(Returns);
rho=(sum(sum(Correlations))-number_stock)/(number_stock*(number_stock-1));
vola= sqrt(diag(Sigma));
Correlations_CC=ones(size(Correlations))*rho+(1-rho)*diag(ones(number_stock,1));
Sigma_CC=diag(vola)*Correlations_CC*diag(vola);
%% 2.4 cov matrix with SHRINKAGE toward COSTANT CORRELATION APPROACH(SCC)
k=0.2;
Sigma_SCC=(1-k)*Sigma*k+k*Sigma_CC;

%% how are located assets on the frontier?
for i=1:number_stock
    hold on
   scatter(sqrt(Sigma_SCC(i,i)),mu_exp(i),"filled")
end

%% 2.6 computing Mark SCC e Mu_exp
A= (ones(number_stock,1)')*(Sigma_SCC\mu);
B=mu'*(Sigma_SCC\mu);
C=ones(number_stock,1)'*(Sigma_SCC\ones(number_stock,1));
D=B*C-A^2;
g=(B*(Sigma_SCC\ones(number_stock,1))-A*(Sigma_SCC\mu))/D;
h=(C*(Sigma_SCC\mu)-A*(Sigma_SCC\ones(number_stock,1)))/D;
m=linspace(A/C,0.1,100);
Var_w=@(m)C/D*(m-A/C).^2+1/C;
hold on
plot(sqrt(Var_w(m)),m,"red")
hold on
%% 2.7 wMVP
w_min=Sigma_SCC\ones(number_stock,1)/C;
mean_MVP = A/C;
scatter(sqrt(1/C),mean_MVP,100,"red","filled")
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 3
%% 3.1 plotto frontiera con sigma_SCC mu_exp e RISK FREE ASSET
Rf=0.01;
A= (ones(number_stock,1)')*(Sigma_SCC\mu);
B=mu'*(Sigma_SCC\mu);
C=ones(number_stock,1)'*(Sigma_SCC\ones(number_stock,1));
K=B-2*A*Rf+C*Rf^2;
m_r=linspace(Rf,0.1,100);
Var_Rf=(m_r-Rf).^2/K;
hold on
plot(sqrt(Var_Rf),m_r)

%% 3.2 Tangent portfolio
P_tang=(Sigma_SCC\(mu_exp-Rf))/sum(Sigma_SCC\(mu_exp-Rf));
%mu_tang
mu_ptang = (A*Rf-B)/(C*Rf-A);
Var_ptang = mu_ptang^2/K;
DV_ptang = sqrt(Var_ptang)
hold on
scatter(DV_ptang,mu_ptang,100,"green","filled")

%% Empiric variance
for i=1:length(m_r)
w=(m_r(i)-Rf)/K*(Sigma_SCC\(mu-Rf));
Var_Rf_emp(i)=w'*Sigma_SCC*w;
end
plot(sqrt(Var_Rf_emp),m_r,'*')
% Numerical solution
Aeq=[mu_exp-Rf]';
for i=1:length(m)
  beq=[m_r(i)-Rf];
  w=(m_r(i)-Rf)/K*(Sigma_SCC\(mu-Rf));
  ww=fmincon(@(ww) ww'*Sigma_SCC*ww,w_min,[],[],Aeq,beq);
  Var_w_num_Rf(i)=ww'*Sigma_SCC*ww;
end
  plot(sqrt(Var_w_num_Rf),m_r,'o');

%%%%%%%%%%%%

hold on
legend("FPstandard","ANF","CATY","FLR","LQDT","NUVA","PII","FP-SCC&mu-exp","wMVP","Var-rf-analitica","P-tang analitico","Var-rf-empirica","Var-rf-numerica")

%% 4.    
%% 4.1.1 Numerical solution no short constraint
%% without rf
% Vincoli
Aeq = [mu_exp'; ones(1, number_stock)]; % Vincolo senza risk-free
% Qui imponiamo ancora la somma ad 1 dei pesi.
lb = zeros(number_stock,1); % Lower bound
ub = ones(number_stock,1); % Upper bound

Var_w_num = zeros(100,1);
m_r=linspace(Rf,0.1,100);

for i = 1:length(m_r)
    beq = [m_r(i), 1]; % somma dei pesi uguale ad 1
    ww = fmincon(@(x) ww' * Sigma_SCC * ww, w_min, [], [], Aeq, beq, lb, ub);
    Var_w_num(i) = ww' * Sigma_SCC * ww;
end
%% With rf asset
Rf=0.01;
A= (ones(number_stock,1)')*(Sigma_SCC\mu_exp);
B=mu_exp'*(Sigma_SCC\mu_exp);
C=ones(number_stock,1)'*(Sigma_SCC\ones(number_stock,1));

K=B-2*A*Rf+C*Rf^2;
m_r=linspace(Rf,0.1,100);
Var_Rf=(m_r-Rf).^2/K; 
figure
hold on
plot(sqrt(Var_Rf),m_r)

Var_w_emp=zeros(size(m));
for i=1:length(m_r)
    w=(m_r(i)-Rf)/K*(Sigma_SCC\(mu_exp-Rf));
    Var_Rf_emp(i)=w'*Sigma_SCC*w;
end
    plot(sqrt(Var_Rf_emp),m_r,'*')
% Numerical solution
Aeq=[mu_exp-Rf]';

for i=1:length(m)
    beq=[m_r(i)-Rf];
    w=(m_r(i)-Rf)/K*(Sigma_SCC\(mu_exp-Rf));
    ww=fmincon(@(ww) ww'*Sigma_SCC*ww,w_min,[],[],Aeq,beq,lb,ub);
    Var_w_num_Rf(i)=ww'*Sigma_SCC*ww;
end
    plot(sqrt(Var_w_num_Rf),m_r,'o');

%% 4.2.2 Difference between return with and without the short-selling constraint
%% With rf asset
Aeq=[mu_exp-Rf]';
lb = [-1;-1;-1;-1;-1;-1];
ub = zeros(number_stock,1);
for i=1:length(m)
    beq=[m_r(i)-Rf];
    w=(m_r(i)-Rf)/K*(Sigma_SCC\(mu_exp-Rf));
    ww=fmincon(@(ww) ww'*Sigma_SCC*ww,w_min,[],[],Aeq,beq,lb,ub);
    Var_w_num_Rf(i)=ww'*Sigma_SCC*ww;
end
    plot(sqrt(Var_w_num_Rf),m_r,'o');

for i=1:length(m)
    beq=[m_r(i)-Rf];
    w=(m_r(i)-Rf)/K*(Sigma_SCC\(mu_exp-Rf));
    ww=fmincon(@(ww) ww'*Sigma_SCC*ww,w_min,[],[],Aeq,beq);
    Var_w_num_Rf(i)=ww'*Sigma_SCC*ww;
end
    plot(sqrt(Var_w_num_Rf),m_r,'o');
%% 4.extra
% % Solving the problem with financial toolbox
% % No short selling constraints
% p1 = Portfolio('AssetMean', mu_exp, 'AssetCovar', Sigma_SCC);
% %% Difference between return with and without the short-selling constraint
% mu=0.005;%expected return
% p2 = Portfolio('AssetMean', mu, 'AssetCovar', Sigma_SCC);
% PortWts2 = estimateFrontier(p2, 1);
% [PortRisk2, PortReturn2] = estimatePortMoments(p2, PortWts2);
% PortRisk2;
% PortReturn2;
% PortWts2;
% % Portfolio with short selling constraints
% AssetBounds = [0,[]]; %in order to deny short sellin lower bound=0
% p3 = Portfolio('AssetMean', mu, 'AssetCovar', Sigma_SCC);
% p3 = setBounds(p3, AssetBounds(1));
% PortWts3 = estimateFrontier(p3, 1);
% [PortRisk3, PortReturn3] = estimatePortMoments(p3, PortWts3);
% PortRisk3;
% PortReturn3;
% PortWts3;

%% 5.
%% 5.1 Market Returns
close = NDXprice(:,"Close");
close = table2array(close);
Market_ret=zeros(size(close)-1);
for i=1:length(close)-1
   Market_ret(i)=(close(i+1)-close(i))./close(i);
end
mean_ret=mean(Market_ret);
var(Market_ret);
used_dates =NDXprice(1:156,"ExchangeDate");
used_dates = table2array(used_dates);
used_dates = datetime(used_dates, "Inputformat", "dd-MMMM-yyyy");
%used_dates
figure
plot(used_dates,Market_ret,"*--","LineWidth",2)
hold on
plot(used_dates,zeros(size(used_dates)),"--",LineWidth=0.5);
hold on
plot(used_dates,ones(size(used_dates)).*mean_ret, "green", LineWidth=2)
%% 5.2 Estimate CAPM alpha and beta
%% 5.2 Compute betas(no alpha)
Mret=transpose(Market_ret(1:155));
betas= zeros(number_stock,1);
for i=1:number_stock
       lm=fitlm(Mret-1,Returns(:,i)-1,"Intercept",false);
       betas(i)=lm.Coefficients.Estimate(1);
end
betas
figure
hist(betas)
xlabel("beta")
ylabel("cardinalità")
%% 5.2.1 Compute betas & alpha
betas= zeros(number_stock,1);
alphas=zeros(number_stock,1);
for i=1:number_stock
       lm=fitlm(Mret-1,Returns(:,i)-1);
       alphas(i)=lm.Coefficients.Estimate(1);
       betas(i)=lm.Coefficients.Estimate(2);
end
betas
alphas
figure
scatter(alphas,betas,"filled")
xlabel("alpha")
ylabel("beta")
%% 6.1 Market Portfolio
A=ANF(:,4);
A=table2array(A);
B=CATY(:,4);
B=table2array(B);
C=FLR(:,4);
C=table2array(C);
D=LQDT(:,4);
D=table2array(D);
E=NUVA(:,4);
E=table2array(E);
F=PII(:,4);
F=table2array(F);
returns=[A,B,C,D,E,F];
%lets calculate MKTCAP assuming that the market is composed by 6 stocks
AA=ANF(155,5);
AA=table2array(AA);
BB=CATY(155,5);
BB=table2array(BB);
CC=FLR(155,5);
CC=table2array(CC);
DD=LQDT(155,5);
DD=table2array(DD);
EE=NUVA(155,5);
EE=table2array(EE);
FF=PII(155,5);
FF=table2array(FF);
market_cap=[AA,BB,CC,DD,EE,FF];
Total_Cap=sum(market_cap);
%we now calculate the weight of the single stock in the market
w_mkt=(market_cap)'/Total_Cap;
%% 6.2 Market Implied Returns π 
% assuming Tau=1
Tau=1;
PI=Tau*sigma*w_mkt;
%here we calculate mean returns to compare them with PI
mean_returns=(mean(returns))';
%% 7
Q=[0.025,0.05]';
certainty=[0.8,0.6];
P=[0,0,0,0,0,1;
   0,-1,1,0,0,0];
%omega matrix 
omega=zeros(size(P,1),size(P,1));
for i=1:size(P,1)
omega(i,i)=P(i,:)*sigma*P(i,:)'*(1-certainty(i))/certainty(i);
end
%here the new sigma
sigma_new=inv(inv(sigma)+P'*inv(omega)*P);
%here the new PI
PI_new= sigma_new*(sigma\PI+P'*inv(omega)*Q);
%lets see the new tangent portfolio
w_new=(sigma_new\PI_new)/sum((sigma_new\PI_new));
%sigma
