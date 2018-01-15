function ImageCreate(disp,prnt)
if nargin < 1, disp = [1,1,1,1,1,1,1,1]; end
if nargin < 2, prnt = [0,0,0,0,0,0,0,0]; end

if disp(1) == 1, FourierTrans(prnt(1));   end
if disp(2) == 1, Fourier_Pulses(prnt(2)); end
if disp(3) == 1, Voronoi_Pol(prnt(3));    end
if disp(4) == 1, Equalization(prnt(4));   end
if disp(5) == 1, Hough(prnt(5));          end
if disp(6) == 1, Edge_Detect(prnt(6));    end
if disp(7) == 1, Noise_Remove(prnt(7));   end
if disp(8) == 1, Degrad(prnt(8));         end
end

% ----------------------------------------------------------------------- %
%% Fourierova transformace
function FourierTrans(prnt)

I{1} = [];                                          % Pruhyvertikalne
    for k=1:5
        I{1} = [I{1},zeros(100,10),ones(100,10)];
    end
I{2} = I{1}';                                       % Pruhy horizontalne
I{3} = [];                                          % Pruhy diagonalne
    pom = [];
    for k=1:5
        pom = [pom,zeros(1,10),ones(1,10)];
    end
    for k=1:100
        I{3} = [I{3};pom];
        pom = [pom(2:end),pom(1)];
    end
I{4} = rot90(I{3});                                 % Pruhy diagonalne
I{5} = [zeros(100,49),ones(100,1),zeros(100,50)];   % Jeden pruh
I{6} = eye(100);                                    % Digonala
I{7} = zeros(100);
I{7}(40:60,45:55) = 1;                              % Obdelnikovy puls
I{8} = kruh(8,100);                                 % Kruhovy puls
I{9} = kruhgauss(5,100);                            % Gaussuv puls

% Vykresleni jednotlivych signalu a jejich amplitudoveho spektra
for i = 1:length(I)
    pom = fftshift(log(abs(fft2(I{i}))+1));
    pom = (pom - min(pom(:)))/max(pom(:));
    
    figure;
    subplot(1,2,1)
    zobr(I{i});
    subplot(1,2,2)
    zobr(pom);
    
    if prnt
        name1 = ['FT_vzor' num2str(i) '.png'];
        name2 = ['FT_obraz' num2str(i) '.png'];
        imwrite(I{i},name1)
        imwrite(pom,name2,'BitDepth',4)
    end
end
end

% Kruh
function K = kruh (R, N)
% vraci binarni kruhovou masku o polomeru R v matici NxN
[X,Y] = meshgrid(-(N-1)/2:(N-1)/2, -(N-1)/2:(N-1)/2);
K = double(X.^2 + Y.^2 < R^2);
end

% Kruh Gauss
function K = kruhgauss(R,N)
% vraci gaussovskou masku o rozptylu R v mat. NxN
[ X, Y] = meshgrid(-(N-1)/2:(N-1)/2,-(N-1)/2:(N-1)/2);
K = exp(-0.5*(X.^2+Y.^2)/R.^2);
end

% ----------------------------------------------------------------------- %
%% 3D znazorneni dulezitych pulsu
function Fourier_Pulses(prnt)
I{1} = zeros(100);
I{1}(40:60,45:55) = 1;                              % Obdelnikovy puls
I{2} = kruh(8,100);                                 % Kruhovy puls
I{3} = kruhgauss(5,100);                            % Gaussuv puls
for k = 1:length(I)
    [surf{k}, org{k}] = ImagePlot3D(I{k}/255);
      
    if prnt
        set(org{k}, 'PaperUnits', 'centimeters' ,'PaperPosition',[0,0,20,20])
        set(surf{k}, 'PaperUnits', 'centimeters','PaperPosition',[0,0,20,20])
        name_1 = ['FT_vzor_3D' num2str(k) '.png'];
        name_2 = ['FT_obraz_3D' num2str(k) '.png'];
        print(org{k},name_1,'-dpng')
        print(surf{k},name_2,'-dpng')
    end
end
end

% 3D plot
function [srf, org] = ImagePlot3D(I)
I = im2double(I);

s = round(size(I)/2);
a = 3;
x = -s(1)+a:s(1)-a;
y = -s(2)+a:s(2)-a;
len_x = length(x);
len_y = length(y);

PS = log(fftshift(abs(fft2(I)))+1);
PS = PS(a:a+len_x-1,a:a+len_y-1);
C = [0,175,0; 0,190,0; 0,205,0; 0,220,0; 0,235,0; 0,250,0]/255;

[Y,X] = meshgrid(y,x);

org = figure;
surf(X,Y,I(a:a+len_x-1,a:a+len_y-1))
axis tight off
colormap(C)

srf = figure;
surf(X,Y,PS);
axis tight off
colormap(C)
end

% ----------------------------------------------------------------------- %
%% Voronoiovy polygony
function Voronoi_Pol(prnt)
x = gallery('uniformdata',[1 20],0);
y = gallery('uniformdata',[1 20],1);
vor = figure;
voronoi(x,y)

if prnt
    set(vor, 'PaperUnits', 'centimeters' ,'PaperPosition',[0,0,20,20])
    print(vor,'VoronoiPolyg','-dpng')
end
end

% ----------------------------------------------------------------------- %
%% Ekvalizovany histogram
function Equalization(prnt)
I{10} = double(imread('Lena_1_Gray.pgm'));
I{11} = linhist(I{10});
I{12} = ekvHist(I{10});

for k = 10:12
    name_5 = ['Lena_obr' num2str(k) '.png'];
    figure;
    zobr(I{k})
    
    ekv = figure;
    histogram(I{k},'Normalization','pdf','FaceColor',[16,158,0]/255)
    xlim([0,255]);
    ylim([0,0.009]);
    xticks([0,255])
    yticks([])
    
    ekv_cum = figure;
    h = histogram(I{k},'Normalization','cdf','FaceColor',[16,158,0]/255);
    xlim([0,255]);
    xticks([0,255])
    yticks([])

    if prnt
        imwrite(I{k}/255,name_5)
        set(ekv, 'PaperUnits', 'centimeters' ,'PaperPosition',[0,0,20,10])
        set(ekv_cum, 'PaperUnits', 'centimeters','PaperPosition',[0,0,20,10])
        name_6 = ['Lena_hist' num2str(k)];
        name_7 = ['Lena_hist_cum' num2str(k)];
        print(ekv,name_6,'-dsvg')
        print(ekv_cum,name_7,'-dsvg')
    end
end
end

% Ekvalizace histogramu
function R = ekvHist(I)
B = 255; % intenzita bile barvy
Vel = length(I(:)); % pocet pixelu v obrazku
R = I; % vysledny snimek
S = 0; % pocet zpracovanych pixelu
for K = 0 : B % prochazime pixely pres intenzity
W = (I == K); % pixely s intenzitou K
P = sum (W(:)); % pocet techto pixelu
R(W) = round((S + P/2) * B / Vel); % nova intenzita
S = S + P; % pocet zpracovanych pixelu
end
end

% Linearni roztazeni histogramu
function R = linhist(I)
% I...obrazek
R = I-min(I(:));
R = R/max(R(:))*255;
end

% ----------------------------------------------------------------------- %
%% Houghova transformace
function Hough(prnt)
I{1} = imnoise(zeros(100),'salt & pepper',0.002);
I{2} = eye(100);

for k = 1:length(I)
    figure;
    zobr(I{k})
    
    [H,~,~] = hough(I{k},'RhoResolution',0.5,'ThetaResolution',0.5);
    figure;
    Hough = imadjust(mat2gray(H));
    imshow(Hough);
    axis off tight
    
    if prnt
        name_1 = ['Hough_org' num2str(k) '.png'];
        name_2 = ['Hough_tra' num2str(k) '.png'];
        imwrite(I{k},name_1)
        imwrite(Hough,name_2)
    end
end
end

%% Detekce hran
function Edge_Detect(prnt)
I = double(imread('Lena_1_Gray.pgm'))/255;
filtr = {'Original','Roberts','Prewitt','Sobel','Canny','Log'};
pom = [1,2,3,4,6,7,8];

figure
for k = 1:length(filtr)
    if pom(k) == 1
        I_f = I;
    else
        I_f = edge(I,filtr{k});
    end
    subplot(2,4,pom(k))
    zobr(I_f)
    title(filtr{k})
    
    if prnt
        name = [filtr{k} '.png'];
        imwrite(I_f,name)
    end
end
end

% ----------------------------------------------------------------------- %
%% Odstraneni sumu
function Noise_Remove(prnt)
I = double(imread('Lena_1_Gray.pgm'))/255;
noise = {imnoise(I,'gaussian',0,0.0025),imnoise(I,'salt & pepper',0.025)};
filtr = {'average','disk','gaussian'};
pom = [1,2,3,4,6,7];

for l = 1:length(noise)
    I_n = noise{l};
    figure
    for k = 1:length(pom)
        if pom(k) == 1
            I_f = I_n;
            ttl = 'Original se sumem';
        elseif pom(k) == 6
            
            I_f = medfilt2(I);
            ttl = 'Median';
        elseif pom(k) == 7
            I_f = wiener2(I_n,[5 5]);
            ttl = 'Wiener';
        else
            I_f = filter2(fspecial(filtr{k-1},5),I_n)/255;
            ttl = filtr{k-1};
        end

        subplot(2,4,pom(k))
        zobr(I_f)
        title(ttl)
        
        if prnt
            name = [ttl '.png'];
            imwrite(I_f,name)
        end
    end
end
end

% ----------------------------------------------------------------------- %
%% Degradace obrazu
function Degrad(prnt)
ttl = {'Original','Rozmazani pohybem','Defokusace','Turbulence',...
    'Rozmazani pohybem FT','Defokusace FT','Turbulence FT'};
pom = [1,2,3,4,6,7,8];

I{1} = double(imread('Lena_1_Gray.pgm'));
I1 = fspecial('motion',10,0); % Obdelnikovy puls
I2 = fspecial('disk',10); %Kruhovy puls
I3 = fspecial('gaussian',10); % Gaussuv puls

I{2} = imfilter(I{1},I1,'replicate');
I{3}= imfilter(I{1},I2,'replicate');
I{4} = imfilter(I{1},I3,'replicate');
I{5} = fftshift(log(abs(fft2(I{2}))+1));
I{6}= fftshift(log(abs(fft2(I{3}))+1));
I{7} = fftshift(log(abs(fft2(I{4}))+1));

figure
for k = 1:length(I)
    subplot(2,4,pom(k))
    zobr(I{k});
    title(ttl{k})
    if prnt
        name = [ttl{k} '.png'];
        imwrite(I{k},name)
    end
end
end

% ----------------------------------------------------------------------- %
%% Pomocne funkce
function zobr(I)
% Img … matice obrazových dat;

colormap(gray(256));
imagesc(I);
axis image;
end
