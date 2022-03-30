% Optionspreis einer Put Option im BS-Modell
% v_put = put_gbm_cos(100,100,0.1,0.25,0.1,0,16)

function v_put = put_gbm_cos(s,K,T,sigma,r,q,N)

L = 10; 
c1 = (r-q)*T; 
c2 = sigma*T; 
a = c1-L*sqrt(c2); 
b = c1+L*sqrt(c2);
x = log(s/K);  
k = (0:N-1)';

Ak = real(charfunc(k*pi/(b-a),x,sigma,r,q,T).*exp(-1i*k*pi*a/(b-a))); 
Vk_put = 2/(b-a)*K*(-chifcn(a,0,a,b,k)+phifcn(a,0,a,b,k));

v_put = exp(-r*T)*(sum(Ak.*Vk_put)-0.5*Ak(1)*Vk_put(1));

% d1 = (log(s/K)+(r-q+sigma^2/2)*T)/(sigma*sqrt(T)); 
% d2 = d1-sigma*sqrt(T);

% bs_call = exp(-q*T)*s*normcdf(d1)-K*exp(-r*T)*normcdf(d2);
% bs_put  = bs_call - exp(-q*T)*s + exp(-r*T)*K; %put/call Parity

%bs_put  = K*exp(-r*T)*normcdf(-d2)-s.*normcdf(-d1);

% e = abs(bs_put-v_put)/bs_put;
end
% ==============================================================================
function f = charfunc(w,x,sigma,r,q,T)

f = exp(1i*w*(x+(r-q-0.5*sigma^2)*T)-0.5*w.^2*sigma^2*T);

end
% ==============================================================================
function f = chifcn(c,d,a,b,k)

f = 1./(1+(k*pi/(b-a)).^2).*(cos(k*pi*(d-a)/(b-a))*exp(d)-...
cos(k*pi*(c-a)/(b-a))*exp(c)+k*pi/(b-a).*sin(k*pi*(d-a)/(b-a))*exp(d)-...
k*pi/(b-a).*sin(k*pi*(c-a)/(b-a))*exp(c));

end
% ==============================================================================
function f = phifcn(c,d,a,b,k)

f = (sin(k*pi*(d-a)/(b-a))-sin(k*pi*(c-a)/(b-a)))*(b-a)./(k*pi);
f(1) = d-c;

end