%1 [CutFrq1] = sp_SetCutOffFreq(x,Fs,'ub',10); uniform band 10

%2 [CutFrq1] = sp_SetCutOffFreq(x,Fs,'d',10); dyadic filter bank 10

%3 [CutFrq1] = sp_SetCutOffFreq(x,Fs,'nd',10,1.8); dyadic filter bank 10

%4 non-dyadic filter bank, Base 'a>1' band filter, e.g. 2.5 or any other real number
% [0 (Fs/2)/a^N ... (Fs/2)/a^4 (Fs/2)/a^3 (Fs/2)/a^2 (Fs/2)/a^1 Fs/2]

%5 [CutFrq1] = sp_SetCutOffFreq(x,Fs,'ee',10); equal energy band 10 

%6 [CutFrq1] = sp_SetCutOffFreq(x,Fs,'elp',10,1); equal Lp-norm band 10

%7 [CutFrq1] = sp_SetCutOffFreq(x,Fs,'pkh',0.2); MinPeakProminence with
%normalized amplitude

function [CutFrq1] = sp_SetCutOffFreq(x,Fs,varargin)

if isempty(varargin{1}) % deafult setup
   NumberOfBands=10;
   %Dyadic=1;
   CutFrq1(1:NumberOfBands+1)=0;
    for i=1:NumberOfBands
        CutFrq1(i)=(Fs/2^(i));
    end
    CutFrq1=fliplr(CutFrq1);
end

if (nargin > 1 && ischar(varargin{1})) && any(strcmpi(varargin{1},{'U','u','UB','ub'}))      
    %UniforBands=1; % enable/disable (1/0) equal frequency bands   
    NumberOfBands=varargin{2};
    CutFrq1(1:NumberOfBands+1)=0;   
    inc=(Fs/2)/NumberOfBands;
    for i=1:NumberOfBands
        CutFrq1(i)=Fs/2-(i-1)*inc;
        %CutFrq1(i)=(i-1)*inc;
    end
    CutFrq1=fliplr(CutFrq1);
end

if (nargin > 1 && ischar(varargin{1})) && any(strcmpi(varargin{1},{'D','d','DB','db'}))      
       %Dyadic, dyadic filter bank
       %[0 (Fs/2)/2^N ... (Fs/2)/2^4 (Fs/2)/2^3 (Fs/2)/2^2 (Fs/2)/2^1 Fs/2]
       NumberOfBands=varargin{2};
       CutFrq1(1:NumberOfBands+1)=0;
    for i=1:NumberOfBands
        CutFrq1(i)=(Fs/2^(i));
    end
    CutFrq1=fliplr(CutFrq1);
end

if (nargin > 1 && ischar(varargin{1})) && any(strcmpi(varargin{1},{'ND','nd','AB','ab'}))      
       %NonDyadic, non-dyadic filter bank  
       % non-dyadic filter bank, Base 'a>1' band filter, e.g. 2.5 or any other real number
       %[0 (Fs/2)/a^N ... (Fs/2)/a^4 (Fs/2)/a^3 (Fs/2)/a^2 (Fs/2)/a^1 Fs/2]
       NumberOfBands=varargin{2};
       CutFrq1(1:NumberOfBands+1)=0;
       a=varargin{3};
       %CutFrq1(1)=Fs/2;
    for i=1:NumberOfBands
        CutFrq1(i)=(Fs/2)/(a^(i-1));
    end
    CutFrq1=fliplr(CutFrq1);
end

if (nargin > 1 && ischar(varargin{1})) && any(strcmpi(varargin{1},{'EE','ee','EEB','eeb'}))      
       %EqualEnergyBand=1; % enable/disable (1/0) dyadic filter bank 
       NumberOfBands=varargin{2};
       CutFrq1(1:NumberOfBands+1)=0;
       
       %%
       L=length(x);
       r=rem(L,2);
       if (r~=0)
           L=L+1; % make it even
       end
       Y = fft(x,L)/L;
       %Compute the two-sided spectrum P2. Then compute the single-sided spectrum P1 based on P2
       %and the even-valued signal length L.
    P2 = (abs(Y)); 
    P1 = P2(1:L/2+1);
    P1(2:end-1) = 2*P1(2:end-1);
    P1=(P1.^2); % square of L2 norm, i.e., energy
    TotalEnergy=sum(P1);
    OneBandEnergy=TotalEnergy/NumberOfBands;
    f = Fs*(0:(L/2))/L;
    %%
    Ex=0; m=2;   
    for k=1:length(P1)
        Ex=Ex+P1(k);
        if((Ex>OneBandEnergy)|| (k==length(P1)))
            CutFrq1(m)=f(k);
            Enrg(m)=Ex;
            Ex=0; % reinitialise
            m=m+1; 
        end
    end
   %tt=1;  
end
if (nargin > 1 && ischar(varargin{1})) && any(strcmpi(varargin{1},{'ELp','elp','ELpB','elpb'}))      
       NumberOfBands=varargin{2};
       p=varargin{3};
       CutFrq1(1:NumberOfBands+1)=0;
       
       %%
       L=length(x);
       r=rem(L,2);
       if (r~=0)
           L=L+1; % make it even
       end
       Y = fft(x,L)/L;
       %Compute the two-sided spectrum P2. Then compute the single-sided spectrum P1 based on P2
       %and the even-valued signal length L.
    P2 = (abs(Y)); 
    P1 = P2(1:L/2+1);
    P1(2:end-1) = 2*P1(2:end-1);
    P1=(P1.^p); % equal Lp norm, i.e., sum of equal absolute values
    TotalEnergy=sum(P1);
    OneBandEnergy=TotalEnergy/NumberOfBands;
    f = Fs*(0:(L/2))/L;
    %%
    Ex=0; m=2;   
    for k=1:length(P1)
        Ex=Ex+P1(k);
        if((Ex>OneBandEnergy)|| (k==length(P1)))
            CutFrq1(m)=f(k);
            Enrg(m)=Ex;
            Ex=0; % reinitialise
            m=m+1; 
        end
    end
   %tt=1;  
end

if (nargin > 1 && ischar(varargin{1})) && any(strcmpi(varargin{1},{'fft1','pkp','FFT1','PKP'}))      
    MinPkProminence=varargin{2};    
    L=length(x); NFFT=2^nextpow2(L);
    Y = fft(x,NFFT)/L; 
    f = Fs/2*linspace(0,1,NFFT/2+1);  
    ampf=2*abs(Y(1:(NFFT/2)+1)); 
    ampf=ampf/max(ampf); % normalize
    %plot(f,ampf);
    [pks,locs] = findpeaks(ampf,f,'MinPeakProminence',MinPkProminence); % 0.3 is OK
    CutFrq1=(locs(1:end-1)+locs(2:end))/2;
    CutFrq1=[0 CutFrq1 Fs/2];
    %tt=1; 
end

if (nargin > 1 && ischar(varargin{1})) && any(strcmpi(varargin{1},{'fft2','pkh','FFT2','PKH'}))      
    pkHight=varargin{2};
    L=length(x); NFFT=2^nextpow2(L);
    Y = fft(x,NFFT)/L; 
    f = Fs/2*linspace(0,1,NFFT/2+1);  
    ampf=2*abs(Y(1:(NFFT/2)+1)); 
    ampf=ampf/max(ampf); % normalize
    %plot(f,ampf);
    [pks,locs] = findpeaks(ampf,f,'MinPeakHeight',pkHight) % pkHight e.g. 20% of Max hight 1, i.e. 0.2 
    CutFrq1=(locs(1:end-1)+locs(2:end))/2;
    CutFrq1=[0 CutFrq1 Fs/2];
    %tt=1; 
end





