clear
clc
runtests(fullfile('Terzaghi','testTerzaghi.m'));
runtests(fullfile('SubDomains','testSubDomains.m'));
runtests(fullfile('MortarConvergence','testMortarPoisson.m'));
%runtests(fullfile('Richards','testRichards.m'));