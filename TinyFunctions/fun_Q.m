function [Q] =  fun_Q(F_0x,F_a,a2disp,Aeff,Order)
eps_0 = evalin('base','eps_0');

atFN{1}   = 'a(t)'; %note: fun_QDelt does not have this, because we usually call this first regardless
for n  = 1:Order                                
avecFN{n} = ['a',num2str(2*n+1)];
atFN{n+1}   = ['a(t)^',num2str(2*n+1)]; 
end

QSUM = 0;
for n = 1:Order
ATermsQ.(avecFN{n}) = Aeff*eps_0*F_0x.*(F_0x/F_a).^(2*(n-1)+2)*a2disp(n);             
QSUM = QSUM+ATermsQ.(avecFN{n});
end
assignin('base','avecFN',avecFN);
assignin('base','atFN',atFN);
assignin('base','ATermsQ',ATermsQ);
Q = QSUM; clear QSUM;
end