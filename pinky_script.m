theta_ph = 0;
theta_el = pi/2;
mu_ph = 0;
mu_el = 1;
sigma_ph = .005;
sigma_el = .005;
k_ph = 0;
k_el = 0;

mesh_s = .01;
bound = 1;

Xin = -bound:mesh_s:bound;
Yin = -bound:mesh_s:bound;

Xmat = ones(length(Yin),1)*Xin;
Ymat = Yin'*ones(1,length(Xin));
Dist = zeros(2*bound*(1/mesh_s)+1,2*bound*(1/mesh_s)+1);

for i=1:(2*bound*(1/mesh_s)+1)
    for j=1:(2*bound*(1/mesh_s)+1)
        Dist(i,j) = j00(psiMinusMinusR((i-(bound*(1/mesh_s)+1))*mesh_s,(j-(bound*(1/mesh_s)+1))*mesh_s,theta_ph,theta_el,sigma_ph,sigma_el,k_ph,k_el,mu_ph,mu_el),psiMinusPlusR((i-(bound*(1/mesh_s)+1))*mesh_s,(j-(bound*(1/mesh_s)+1))*mesh_s,theta_ph,theta_el,sigma_ph,sigma_el,k_ph,k_el,mu_ph,mu_el),psiPlusMinusR((i-(bound*(1/mesh_s)+1))*mesh_s,(j-(bound*(1/mesh_s)+1))*mesh_s,theta_ph,theta_el,sigma_ph,sigma_el,k_ph,k_el,mu_ph,mu_el),psiPlusPlusR((i-(bound*(1/mesh_s)+1))*mesh_s,(j-(bound*(1/mesh_s)+1))*mesh_s,theta_ph,theta_el,sigma_ph,sigma_el,k_ph,k_el,mu_ph,mu_el));
    end
end

N=10000;
res=10;
vals=zeros(2,N);
for i=1:N
[vals(1,i),vals(2,i)]=pinky(Xin,Yin,Dist,res);
end

figure;
hold on
imagesc(Xin,Yin,Dist)
colormap(gray)
plot(vals(1,:),vals(2,:),'r.')
xlabel('Xin')
ylabel('Yin')
axis equal tight
box on

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
function T = gendist(P,N,M,varargin)
%
%GENDIST - generate random numbers according to a discrete probability
%distribution.
%Tristan Ursell, (c) 2012
%
%T = gendist(P,N,M)
%T = gendist(P,N,M,'plot')
%
%The function gendist(P,N,M) takes in a positive vector P whose values
%form a discrete probability distribution for the indices of P. The
%function outputs an N x M matrix of integers corresponding to the indices
%of P chosen at random from the given underlying distribution.
%
%P will be normalized, if it is not normalized already. Both N and M must
%be greater than or equal to 1.  The optional parameter 'plot' creates a
%plot displaying the input distribution in red and the generated points as
%a blue histogram.
%
%Conceptual EXAMPLE:
%
%If P = [0.2 0.4 0.4] (note sum(P)=1), then T can only take on values of 1,
%2 or 3, corresponding to the possible indices of P.  If one calls 
%gendist(P,1,10), then on average the output T should contain two 1's, four
%2's and four 3's, in accordance with the values of P, and a possible 
%output might look like:
%
%T = gendist(P,1,10)
%T = 2 2 2 3 3 3 1 3 1 3
%
%If, for example, P is a probability distribution for position, with
%positions X = [5 10 15] (does not have to be linearly spaced), then the
%set of generated positions is X(T).
%
%EXAMPLE 1:
%
%P = rand(1,50);
%T = gendist(P,100,1000,'plot');
%
%EXAMPLE 2:
%
%X = -3:0.1:3;
%P = 1+sin(X);
%T = gendist(P,100,1000,'plot');
%Xrand = X(T);
%
%Note:
%size(T) = 100 1000
%size(Xrand) = 100 1000
%
%Thanks to Derek O'Connor for tips on speeding up the code.
%check for input errors
if and(nargin~=3,nargin~=4)
    error('Error:  Invalid number of input arguments.')
end
if min(P)<0
    error('Error:  All elements of first argument, P, must be positive.')
end
if size(P,1)>1
    P=P';
end
if or(N<1,M<1)
    error('Error:  Output matrix dimensions must be greater than or equal to one.')
end
%normalize P
Pnorm=[0 P]/sum(P);
%create cumlative distribution
Pcum=cumsum(Pnorm);
%create random matrix
N=round(N);
M=round(M);
R=rand(1,N*M);
%calculate T output matrix
V=1:length(P);
[~,inds] = histc(R,Pcum); 
T = V(inds);
%shape into output matrix
T=reshape(T,N,M);
%if desired, output plot
if nargin==4
    if strcmp(varargin{1},'plot')
        Pfreq=N*M*P/sum(P);
        figure;
        hold on
        hist(T(T>0),1:length(P))
        plot(Pfreq,'r-o')
        ylabel('Frequency')
        xlabel('P-vector Index')
        axis tight
        box on
    end
end
end