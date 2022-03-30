% Optionspreis einer Call Option im Heston-Modell mit COS-Methode
% vc_heston = call_heston_cos(100,100,0.5751,0,0,1.5768,0.0398,-0.5711,0.0175,1,16)
% andere parameter:
% s0 = 100; K = 100; r = 0.03; T = 0.5; theta = 0.04;
% v0 = 0.04; q = 0; rho = -0.7; kappa = 2; sigma = 0.5;

function vc_heston = call_heston_cos(s0,K,sigma,r,q,kappa,theta,rho,v0,T,N)

L = 12;
c1 = (r-q)*T+(1-exp(-kappa*T))*((theta-v0)/2*kappa)-0.5*theta*T; 
c2 = 1/(8*kappa^3)*(sigma*T*kappa*exp(-kappa*T)*(v0-theta)*(8*kappa*rho-4*sigma)...
     +kappa*rho*sigma*(1-exp(-kappa*T))*(16*theta-8*v0)...
     +2*theta*kappa*T*(-4*kappa*rho*sigma+sigma^2+4*kappa^2)...
     +sigma^2*((theta-2*v0)*exp(-2*kappa*T)+theta*(6*exp(-kappa*T)-7)+2*v0)...
     +8*kappa^2*(v0-theta)*(1-exp(-kappa*T)));
a = c1-L*sqrt(abs(c2));
b = c1+L*sqrt(abs(c2));
x = log(s0/K); %für mehrere K, schreibe log(s0./K)
k = (0:N-1)';

Ak = real(charfunc(k*pi/(b-a),x,sigma,r,q,kappa,theta,rho,v0,T).*...
     exp(-1i*k.*pi*a/(b-a)));
Vk = 2/(b-a)*K*(chifcn(0,b,a,b,k)-phifcn(0,b,a,b,k));

vc_heston = exp(-r*T)*(sum(Ak.*Vk)-0.5*Ak(1)*Vk(1)); format long;

% ========== Call Preis durch Put-Call Parität: ==========
%vp_heston = put_heston_cos(s0,K,sigma,r,q,kappa,theta,rho,v0,T,N);
%vc_heston = (vp_heston)+s0*exp(-q*T)-K*exp(-r*T);format long

%e = abs(5.785155450-vc_heston)
end
% ========================================================================
function f = charfunc(u,x,sigma,r,q,kappa,theta,rho,v0,T)

d = sqrt((rho*sigma*1i*u-kappa).^2+sigma^2*(1i*u+u.^2));
g2 = (kappa-rho*sigma*1i*u-d)./(kappa-rho*sigma*1i*u+d);

f = exp(1i*u*(x+(r-q)*T));
f = f.*exp(theta*kappa/sigma^2*((kappa-rho*sigma*1i*u-d)*T-...
    2*log((1-g2.*exp(-d*T))./(1-g2))));
f = f.*exp(v0/sigma^2*(kappa-rho*sigma*1i*u-d).*(1-exp(-d*T))./(1-g2.*exp(-d*T)));

% ============== ACHTUNG HESTON TRAP ===============
% gamma = sqrt(sigma^2*(u.^2 + 1i*u)+(kappa - 1i*rho*sigma*u).^2);
% 
% Zaehler1 = (exp(1i*u*(x+(r-q)*T)+kappa*theta*T*(kappa-1i*rho*sigma*u)/sigma^2));
% Nenner1  = (cosh(gamma*T/2)+(kappa-1i*rho*sigma*u)./gamma.*sinh(gamma*T/2)).^...
%            (2*kappa*theta/sigma^2);
% Zaehler2 = (-(u.^2+1i*u)*v0);
% Nenner2  = gamma.*coth(gamma*T/2)+kappa-1i*rho*sigma*u;
% 
% Element1 = (Zaehler1./Nenner1);
% Element2 = (exp(Zaehler2./Nenner2));
% 
% f = Element2.*Element1;
% ============== ACHTUNG HESTON TRAP ===============
end
% ========================================================================
function f = chifcn(c,d,a,b,k)

f = 1./(1+(k*pi/(b-a)).^2).*(cos(k*pi*(d-a)/(b-a))*exp(d)-...
cos(k*pi*(c-a)/(b-a))*exp(c)+k*pi/(b-a).*sin(k*pi*(d-a)/(b-a))*exp(d)-...
k*pi/(b-a).*sin(k*pi*(c-a)/(b-a))*exp(c));
  
end
% ========================================================================
function f = phifcn(c,d,a,b,k)

f = (b-a)./(pi*k).*(sin(k*pi*(d-a)/(b-a))-sin(k*pi*(c-a)/(b-a)));
f(1) = d-c;

end