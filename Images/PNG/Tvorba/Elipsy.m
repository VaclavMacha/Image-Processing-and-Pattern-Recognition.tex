clear; close all; clc

addpath('C:\Users\Vaclav\OneDrive\Vaclav\Dokumenty\FJFI\5. rocnik\ROZ\Cvika\Functions')


% % Primka
% mu_1 = [0 0];
% sig_1 = [0.25 0.3; 0.3 1];
% mu_2 = [2 0];
% sig_2 = [0.25 0.3; 0.3 1];

% % Kruznice
% mu_1 = [0 0];
% sig_1 = [4 0; 0 4];
% mu_2 = [0 0];
% sig_2 = [0.1 0; 0 0.1];

% % Hyperbola
% mu_1 = [0 0];
% sig_1 = [1 0; 0 4];
% mu_2 = [-5 0];
% sig_2 = [4 0; 0 1];

% % Rovnobezky
% mu_1 = [0 0];
% sig_1 = [1 0; 0 1];
% mu_2 = [-4 0];
% sig_2 = [4 0; 0 1];

% % Elipsa
% mu_1 = [0 0];
% sig_1 = [0.4 0; 0 0.1];
% mu_2 = [-4 0];
% sig_2 = [4 0; 0 1];

% Kruznice
mu_1 = [0 0];
sig_1 = [0.1 0; 0 0.1];
mu_2 = [-4 0];
sig_2 = [4 0; 0 1];


x_min = -10;
x_max = 10;

y_min = -10;
y_max = 10;

x = linspace(x_min, x_max);
y = linspace(y_min, y_max);

[X,Y] = meshgrid(x,y);
s = size(X);
G_1 = mvnpdf([X(:) Y(:)],mu_1,sig_1);
G_1 = reshape(G_1,s);

G_2 = mvnpdf([X(:) Y(:)],mu_2,sig_2);
G_2 = reshape(G_2,s);

m1 = mu_1';
m2 = mu_2';
f = @(x,y) ([x;y] - m1)'*inv(sig_1)*([x;y] - m1) + log(det(sig_1))...
    - ([x;y] - m2)'*inv(sig_2)*([x;y] - m2) - log(det(sig_1));

C1 = [0,175,0; 0,190,0; 0,205,0; 0,220,0; 0,235,0; 0,250,0]/255;
C2 = [175,0,0; 190,0,0; 205,0,0; 220,0,0; 235,0,0; 250,0,0]/255;


figure;

hold on
h1 = surf(x,y,G_1);
[~,hc1] = contour(x,y,G_1,1);
hc1.ContourZLevel = -1;
colormap(C1);
freezeColors

surf(x,y,G_2)
[~,hc2] = contour(x,y,G_2,1);
hc2.ContourZLevel = -1;
colormap(C2);

h = ezplot(f);
set(h, 'LineStyle', '--','Color','b')

hold off
axis tight
grid on
view(3)




















