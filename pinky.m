function [x0,y0]=pinky(Xin,Yin,dist_in,varargin)
%check input
if length(size(dist_in))>2
    error('The input must be a N x M matrix.')
end
%check sizes
[sy,sx]=size(dist_in);
if or(length(Xin)~=sx,length(Yin)~=sy)
    error('Dimensions of input vectors and input matrix must match.')
end
%check values
if any(dist_in(:)<0)
    error('All input probability values must be positive.')
end
%get res
if nargin==4
    res=varargin{1};
    if res<=1
        error('The resolution factor (res) must be an integer greater than one.')
    end
elseif nargin~=3
    error('Incorrect number of input arguments.')
end
    
%create column distribution and pick random number
col_dist=sum(dist_in,1);
%pick column distribution type
if nargin==3;
    %if no res parameter, simply update X/Yin2
    col_dist=col_dist/sum(col_dist);
    Xin2=Xin;
    Yin2=Yin;
else
    %generate new, higher res input vectors
    Xin2=linspace(min(Xin),max(Xin),round(res*length(Xin)));
    Yin2=linspace(min(Yin),max(Yin),round(res*length(Yin)));
    
    %generate interpolated column-sum distribution
    col_dist=interp1(Xin,col_dist,Xin2,'pchip');
    
    %check to make sure interpolated values are positive
    if any(col_dist<0)
        col_dist=abs(col_dist);
        warning('Interpolation generated negative probability values.')
    end
    col_dist=col_dist/sum(col_dist);
end
%generate random value index
ind1=gendist(col_dist,1,1);
%save first value
x0=Xin2(ind1);
%find corresponding indices and weights in the other dimension
[val_temp,ind_temp]=sort((x0-Xin).^2);
if val_temp(1)<eps %if we land on an original value
    row_dist=dist_in(:,ind_temp(1));
else %if we land inbetween, perform linear interpolation
    low_val=min(ind_temp(1:2));
    high_val=max(ind_temp(1:2));
    
    Xlow=Xin(low_val);
    Xhigh=Xin(high_val);
    
    w1=1-(x0-Xlow)/(Xhigh-Xlow);
    w2=1-(Xhigh-x0)/(Xhigh-Xlow);
    
    row_dist=w1*dist_in(:,low_val) + w2*dist_in(:,high_val);
end
%pick column distribution type
if nargin==3;
    row_dist=row_dist/sum(row_dist);
else
    row_dist=interp1(Yin,row_dist,Yin2,'pchip');
    row_dist=row_dist/sum(row_dist);
end
%generate random value index
ind2=gendist(row_dist,1,1);
%save first value
y0=Yin2(ind2);
end