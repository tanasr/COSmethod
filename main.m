
% ========== Europäische Call im Black Scholes Modell ==========
s0=100; %Aktueller Preis Underlying
K=100; %Ausübungspreis (strike)
T=0.1; % Zeit bis Ausübung
r=0.1; %Risfree rate
q=0; % Dividende

% Schätzung
sigma=0.25; %Vola

% für COS method
N=128;



vcall_cos = call_gbm_cos(s0,K,T,sigma,r,q,N); %value cos method
vput_cos = put_gbm_cos(s0,K,T,sigma,r,q,N);
[vcall_bs, vput_bs] = blsprice(s0,K,T,r,sigma); %value black scholes
err = sqrt((vcall_bs - vcall_cos)^2);

disp('Europäische Call im Black Scholes Modell')
fprintf('Cos methode : %.9f\n', vcall_cos)
fprintf('Analytisch  : %.9f\n', vcall_bs)
fprintf('Error mit N = %i\n', N);
fprintf('Error = %.9f\n', err);

err = sqrt((vput_bs - vput_cos)^2);

disp('Europäische Put im Black Scholes Modell')
fprintf('Cos methode : %.9f\n', vput_cos)
fprintf('Analytisch  : %.9f\n', vput_bs)
fprintf('Error mit N = %i\n', N);
fprintf('Error = %.9f\n', err);
disp('--------------------')




% ========== Discretely monitored Barrier Option im Black-Scholes Modell ==========
% TODO: DEBUG, error in idx for function Mfcn ...


% ========== Europäische Call im Heston Modell ==========
% vc_heston = call_heston_cos(100,100,0.5751,0,0,1.5768,0.0398,−0.5711,0.0175,1,16)
s0 = 100;
K = 100;
sigma = 0.5751;
r = 0;
q = 0;
kappa = 1.5768;
theta = 0.0398;
rho = -0.5711;
v0 = 0.0175;
T = 1;
N = 128;

vc_heston = call_heston_cos(s0,K,sigma,r,q,kappa,theta,rho,v0,T,N);
true = 5.785155450;
err = abs(5.785155450-vc_heston);

disp('Europäische Call im Heston Modell')
fprintf('Cos methode : %.9f\n', vc_heston)
fprintf('Analytisch  : %.9f\n', true)
fprintf('Error mit N = %i\n', N);
fprintf('Error = %.9f\n', err);
disp('--------------------')



% some tests
KK = (20:2:150)';
n = length(KK);
for j = 1:n
   v(j,1)=call_gbm_cos(100,KK(j,1),0.1,0.25,0.1,0,32);
end

plot(KK,v)
title(['Preis Call Option mit s0 = ',num2str(s0)]);
ylabel('Preis');
xlabel('Strike');






