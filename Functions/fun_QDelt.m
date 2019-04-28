function Q = fun_QDelt(F_0x,F_0y,F_a,a2disp,Aeff,Order)
eps_0 = evalin('base','eps_0');

%                                  eps_0.*F_0x.*(F_0y/F_a).^2*(1/3 * a2n(1)+...
%                                                              1/5 * a2n(2).*(F_0y./F_a).^2+...
%                                                              1/7 * a2n(3).*(F_0y./F_a).^4+...
%                                                              1/9 * a2n(4).*(F_0y./F_a).^6+...
%                                                              1/11* a2n(5).*(F_0y./F_a).^8)*Aeff;
                                                         
                                                         

for n  = 1:Order                                
avecFN{n} = ['a',num2str(2*n+1)];
end

QSUM = 0;
for n = 1:Order
ATermsQDelt.(avecFN{n}) = (1/(2*n+1))*Aeff*eps_0*F_0x.*(F_0y/F_a).^(2*(n-1)+2)*a2disp(n);             
QSUM = QSUM+ATermsQDelt.(avecFN{n});
end
assignin('base','avecFN',avecFN);
assignin('base','ATermsQDelt',ATermsQDelt);
Q = QSUM; clear QSUM;
end