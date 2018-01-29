function Focal Adhesion_analysis()
clear; clc;
originaldir = pwd;
% cd('/Users/admin/Desktop/');
cd('/Volumes/JHC External/JHC Lodoen Data/2017.8/8.7.17 - vinculin coloc');
% cd('/Volumes/JHC External/JHC Lodoen Data/2017.7/7.21.17 - pax coloc');
[FileName,PathName,FilterIndex] = uigetfile({'*.tif'},'MultiSelect','on');
cd(PathName);
CSVname = getNameGUI();
exportCSV(CSVname, [{'Cell'},{'Pearsons'},{'Manders'},{'B1b CTCF'},{'FAK CTCF'},{'Area'}]);
[~,nFiles] = size(FileName);
fileCounter = 0;
for img=FileName
    fileCounter = fileCounter + 1;
    disp("File " + num2str(fileCounter) + " of " + num2str(nFiles));
    disp(img);
    file = imread(char(img));
    [B1backgroundMFI,FAKbackgroundMFI] = background(file); % get background MFI for CTCF
    disp("File " + num2str(fileCounter) + " of " + num2str(nFiles));
    disp(img);
    figure, imshow(file); % for free-form drawing
    contin = 1;
    counter = 0;
    while contin
        position = [];
        while isempty(position)
            form = imfreehand; % outline the cell
            position = wait(form); % double click on shape to continue
        end
        
        mask(:,:,1) = createMask(form);
        mask(:,:,2) = createMask(form);
        mask(:,:,3) = createMask(form);
        cell = immultiply(file, mask);
        
        % send for further calculations
        [pearson, manders] = colocalizationCalcs(cell, mask(:,:,1));
        [~, B1intCTCF] = correctedTotalCellularFluorsence(cell(:,:,1), B1backgroundMFI, mask(:,:,1));
        [Area, FAKCTCF] = correctedTotalCellularFluorsence(cell(:,:,2), FAKbackgroundMFI, mask(:,:,1));
        contin = questionDialogBox();  % 0 = next image, 1 = same image
        
        % save to <data_collect>
        counter = counter + 1;
        name = strcat(char(img),'_',num2str(counter),'n');
        celldata = [{name} num2str(pearson) num2str(manders) num2str(B1intCTCF) num2str(FAKCTCF) num2str(Area)];
        exportCSV(CSVname,celldata);
        disp("File " + num2str(fileCounter) + " of " + num2str(nFiles));
        disp(img);
    end
    close;
end
cd(pwd);
disp('done');
end

function [rP, rM] = colocalizationCalcs(cellROI, cellmask)
% obtains the major clocolization values
B1int = cellROI(:,:,1);
FAK = cellROI(:,:,2);
length = lengthOfMatrix(B1int);
B1int_avg = averageInROI(B1int, cellmask, length);
FAK_avg = averageInROI(FAK, cellmask, length);

rP = Pearson(FAK, B1int, FAK_avg, B1int_avg, length);
rM = MandersOverlap(FAK, B1int, length);
end

function r = Pearson(R, G, Ravg, Gavg, len)
numer = 0; % numerator
denomR = 0; % R part of denominator
denomG = 0; % G part of denominator
for i = 1:len
    Ri = double(R(i));
    Gi = double(G(i));
    numer = numer + ((Ri - Ravg)*(Gi - Gavg));
    denomR = denomR + ((Ri - Ravg)^2);
    denomG = denomG + ((Gi - Gavg)^2);
end
denom = denomR * denomG;
denom = sqrt(denom);
r = numer/denom;
end

function o = MandersOverlap(R, G, len)
numer = 0;
denomR = 0;
denomG = 0;
for i = 1:len
    Ri = double(R(i));
    Gi = double(G(i));
    numer = numer + (Ri*Gi);
    denomR = denomR + (Ri^2);
    denomG = denomG + (Gi^2);
end
denom = (denomR*denomG)^0.5;
o = numer/denom;
end

function l = lengthOfMatrix(m)
[x, y] = size(m);
l = x*y;
end

function a = averageInROI(color, bw, len)
values = [];
for i = 1:len
    if bw(i) == 1
        values = [values color(i)];
    end
end
a = mean(values);
end

function [area, ctcf] = correctedTotalCellularFluorsence(B1, backMFI, bw)
% measure the CTCF
scale = 1/(5.544^2);    % ratio of um:pixels
area = sum(bw(:)) * scale;
length = lengthOfMatrix(bw);
cellMFI = averageInROI(B1, bw, length);
ctcf = area * (cellMFI - backMFI);
end

function [bMFI,fMFI] = background(name)
imshow(name);
B1 = name(:,:,1); % save B1 integrin image
FAK = name(:,:,2); % save FAK integrin image
boxes = zeros(3,4); 
for i = 1:3
    h = imrect; % draw a rectangle
    p = wait(h); % and double click on it
    boxes(i,:) = getPosition(h);
end
% get the individual background rectangles of B1 integrin
b1 = B1(boxes(1,2):(boxes(1,2)+boxes(1,4)) , boxes(1,1): (boxes(1,1)+boxes(1,3)));
b2 = B1(boxes(2,2):(boxes(2,2)+boxes(2,4)) , boxes(2,1): (boxes(2,1)+boxes(2,3)));
b3 = B1(boxes(3,2):(boxes(3,2)+boxes(3,4)) , boxes(3,1): (boxes(3,1)+boxes(3,3)));
% get the individual background rectangles of B1 integrin
f1 = FAK(boxes(1,2):(boxes(1,2)+boxes(1,4)) , boxes(1,1): (boxes(1,1)+boxes(1,3)));
f2 = FAK(boxes(2,2):(boxes(2,2)+boxes(2,4)) , boxes(2,1): (boxes(2,1)+boxes(2,3)));
f3 = FAK(boxes(3,2):(boxes(3,2)+boxes(3,4)) , boxes(3,1): (boxes(3,1)+boxes(3,3)));
% make fake images with all ones for the computeColoc function
o1 = uint8(ones(size(b1)));
o2 = uint8(ones(size(b2)));
o3 = uint8(ones(size(b3)));
% get the MFIs for all the background sections
B1M_data = [0,0,0];
[~,B1M_data(1),~] = computeColoc(o1, b1); 
[~,B1M_data(2),~] = computeColoc(o2, b2);
[~,B1M_data(3),~] = computeColoc(o3, b3);
FAKM_data = [0,0,0];
[~,FAKM_data(1),~] = computeColoc(o1, f1); 
[~,FAKM_data(2),~] = computeColoc(o2, f2);
[~,FAKM_data(3),~] = computeColoc(o3, f3);
bMFI = mean(B1M_data); % return average for B1 background MFI
fMFI = mean(FAKM_data); % return average for FAK background MFI
close;
end

% just using this for the background function
function [A,m,int_den] = computeColoc(d, s)
% computes area, MFI, and integrated density
mult = immultiply(d, s);
s_total = sum(mult(:));
A = sum(d(:));
m = s_total / A;
int_den = m*A;
end

function r = questionDialogBox()
answer = questdlg('Would you like to continue with this image?',...
    '',...
    'Yes','No','Yes');
switch answer
    case 'Yes'
        r = 1;
    case 'No'
        r = 0;
end
end

function exportCSV(fname, d)
% brute force to make add to CSV file
if isempty(d) ~= 1
    fid = fopen(fname,'a');
    if fid>0
        for k=1:size(d,1)
            fprintf(fid,'%s,%s,%s,%s,%s,%s\n',d{k,:});
        end
    end
    fclose(fid);
end
end

function n = getNameGUI()
% name the CSV file
prompt = {'Name csv file'};
dlg_title = 'Save as';
defaultans = {''};
a = inputdlg(prompt,dlg_title,[1 50],defaultans);
n = a{1};
n = strcat(n,'.csv');
end