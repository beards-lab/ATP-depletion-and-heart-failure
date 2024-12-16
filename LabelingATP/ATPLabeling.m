kL   = 5;
kS1T = 10;
kS2T = 0.2;
kS1D = 0.1;
kS2D = 10;
kH   = 5;
k = [kS1T, kS2T, kH, kS1D, kS2D, kL];
 
% setting up and solving linear steady-state problem
A = [-kS2T, +kS1T,        0, +kL;
     +kS2T, -kH-kS1T,     0,   0;
         0, +kH, -kS1D,   kS2D; % fix typo
        +1,  +1,         +1,  +1];
b = [0; 0; 0; 1];
p = A\b;
 
% The label kinetic problem
f0 = [1 1 1 1];
[t,f] = ode15s(@dfdt1,[0 600],f0,[],p,k);
label = f*p;
figure(1); semilogy(t,label)

function g = dfdt1(~,f,p,k)
 
% p = [pST, pUT, pUD, pSD], the steady-state probs.
% k = [kS2T, kH, kS1D, kL], the rate constants
 
fST = f(1);
fUT = f(2);
fUD = f(3);
fSD = f(4);

pST = p(1);
pUT = p(2);
pUD = p(3);
pSD = p(4);

kS1T = k(1);
kS2T = k(2);
kH   = k(3);
kS1D = k(4);
kS2D = k(5);
kL   = k(6);
 
dfST = (-kS2T*fST*pST + kS1T*fUT*pUT)/pST;
dfUT = (-kH*pUT*fUT - kS1T*fUT*pUT + kS2T*fST*pST)/pUT;
dfUD = (-kS1D*pUD*fUD + kH*fUT*pUT + kS2D*fSD*pSD)/pUD;
dfSD = (-kL*pSD*fSD + kS1D*fUD*pUD - kS2D*fSD*pSD)/pSD;
 
g = [dfST; dfUT; dfUD; dfSD];
end