function im = Imgadjust(im , thresh , ol , or)
if nargin<3
    ol=0.01;
    or = 0.99;
end

TypeList = {'uint8' , 'uint16'};
Type = class(im);
idx_class = strcmpi(Type , TypeList);
if isempty(find(idx_class,1))
    err = 'acceptance data for,ats for the input image are 8-bit and 16-bit';
    error(err)
end
if find(idx_class,1)==1
    coef =1;
    BPP = 8;
else
    coef = 2^8;
    BPP = 16;
end
for ch = 1:size(im,3)
    img = im(:,:,ch);
    skew = skewness(double(img(:)));
    c = imhist(img);
    c1 = imerode(c, strel('line' , 5, 90));
    maxima = find(c1 == imdilate(c1 , strel('line' , 90 , 90))  & c1>size(im,1) *size(im,2)/255);
    maxima = maxima(find(maxima > 5 & maxima <250));
    if length(maxima)>1
        if maxima(2)-maxima(1) >100
            x = [maxima(1) maxima(2)];
            y = [maxima(2) maxima(2)-(maxima(2)-maxima(1))/2];
            p = polyfit(x,y,1);
            if BPP==16
                img = uint16(ppval(p,double(img)));
            else
                img = uint8(ppval(p , double(img)));
            end
        end
    end
    counts = imhist(img);
    counts=medfilt1(counts , 5);
    counts(1) = 0;
    cum_sum = cumsum(counts);
    if skew > 0.6
        threshL = 0.05*thresh;
        threshR = 0.06*thresh;
    elseif skew>0.25 & skew<=0.06
        threshL = 0.1*thresh;
        threshR = 0.6*thresh;
    elseif skew>-0.25 & skew<=0.25
        threshL = 0.1*thresh;
        threshR = 0.01*thresh;
    elseif skew>-0.6 & skew<=-0.25
        threshL = 0.5*thresh;
        threshR = 0.3*thresh;
    else
        threshL = 0.5*thresh;
        threshR = 0.03*thresh;
    end
    Left = find(cum_sum < threshL * cum_sum(end));
    Left = length(Left)/256;
    cum_sum = cumsum(fliplr(counts'))';
    Right = find( cum_sum <threshR * cum_sum(end));
    Right = length(Right)/256;
    a = img<Left *2^BPP;
    b = img> (1-Right) * 2^BPP;
    a1 = double(a) .* double(img);
    b1 = double(b) .* double(img);
    a1 = a1 .* (0.03/Left);
    b2 = (b1 - (1- Right)*2^BPP) *(2^BPP - 0.98 * 2^BPP) / (2^BPP - (1 - Right) * 2^BPP)+ 0.98 * 2^BPP;
    img1 = imadjust(img , [Left 1-Right] , [ol or]);
    img1(a) = a1(a);
    img1(b) = b2(b);
    im(:,:,ch) = img1;
end
end








