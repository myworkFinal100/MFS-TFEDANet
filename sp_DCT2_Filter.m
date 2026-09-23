function FilteredData1=sp_DCT2_Filter(data_xn,Fs,varargin) % FFT based Intrinsic Mode Decomposition 
    FigureEnable=0;
    gaindB=0;
    [RowData,ColomnData] = size(data_xn); % change data to row vector   
    if(RowData>ColomnData)
        data_xn=data_xn';   
    end    
    N=length(data_xn);
    Xk=dct(data_xn); % compute DCT2 
            
    
    if FigureEnable
        freq = (Fs/2)*linspace(0,1,N);
        DataFFT=abs(Xk);        
        subplot(2,1,1);
        if gaindB
        plot(freq,20*log10(DataFFT)); 
        else
        plot(freq,(DataFFT)); 
        end
        title('Data DCT plot');
    end    
    
    %fc=Fs*k/(2*N) => k=fc*N*2/Fs
    %range_k=ceil(fc*N*2/Fs);
    %range_k=round(fc*N*2/Fs); % convert fc to k (digital freq)
    Hk=zeros(1,N);
    if (nargin > 1 && ischar(varargin{end})) && any(strcmpi(varargin{end},{'LP','lp'})),
        if strcmpi(varargin{end},'LP')
            fc=varargin{end-1};
            range_k=round(fc*N*2/Fs); % convert fc to k (digital freq)
            kk{1}=1:range_k; % for low pass filter             
            Hk(kk{1})=1;            
            Xk((range_k+1):N)=0;            
        end        
    end
    
    if (nargin > 1 && ischar(varargin{end})) && any(strcmpi(varargin{end},{'HP','hp'})),
        if strcmpi(varargin{end},'HP')
            fc=varargin{end-1};
            range_k=round(fc*N*2/Fs); % convert fc to k (digital freq)            
            kk{1}=(range_k+1):N; % for High pass filter              
            Hk(kk{1})=1;            
            Xk(1:range_k)=0;            
        end        
    end
    if (nargin > 1 && ischar(varargin{end})) && any(strcmpi(varargin{end},{'BP','bp'})),
        if strcmpi(varargin{end},'BP')
            fc1=varargin{end-2};
            fc2=varargin{end-1};
            range_k1=round(fc1*N*2/Fs); % convert fc to k (digital freq) 
            range_k2=round(fc2*N*2/Fs); % convert fc to k (digital freq)
            kk{1}=(range_k1+1):range_k2; % for High pass filter   
              
            Hk(kk{1})=1;        
            
            Xk(1:range_k1)=0;
            Xk(range_k2+1:N)=0;           
        end         
    end 
    
    if (nargin > 1 && ischar(varargin{end})) && any(strcmpi(varargin{end},{'BS','bs'}))
        if strcmpi(varargin{end},'BS')
            Hk=ones(1,N);
            fc1=varargin{end-2};
            fc2=varargin{end-1};
            range_k1=round(fc1*N*2/Fs); % convert fc to k (digital freq) 
            range_k2=round(fc2*N*2/Fs); % convert fc to k (digital freq)
            kk{1}=(range_k1+1):range_k2; % for High pass filter   
            Hk(kk{1})=0;            
            Xk(range_k1:range_k2)=0;                           
        end        
    end    
    
    if FigureEnable
        subplot(2,1,2);
        plot(freq,abs(Xk));
    end
    
    if 0
        
        if 0 % direct method which computationally complex
            n=1:N;
            w=ones(1,N)*(2/N)^(0.5); w(1)= (1/N)^(0.5); 
            tmp=0;
            for k=1:N        
                    tmp=tmp+w(k)*Xk(k)*exp( 1i*pi*(2*n-1)*(k-1)/(2*N) ); % DCT-2 and DCT-3 are inverses of each other.
            end            
            FilteredData=tmp; % retun analytic signal   
        else % FFT based               
            % IDCT2 using 4N FFT
            ak=(1/sqrt(N)) * [1; sqrt(2)*ones(N-1,1)]; % Normalization factor
            ck=ak'.*Xk;
            fx4N=fft(ck,4*N);           
            xfsas=fx4N(2:2:(2*N));% FSAS 

            xr=real(xfsas);   % signal x         
            xq=-imag(xfsas);  % quadrature component, - sign due to -j in exp(-j*2*pi*k*n/(N))
            FilteredData= xr+1i*xq; % retun Fourier-Singh analytic signal (FSAS)           
        end
        
    else
        %IMFs_y=ifft(Hk.*Xk)*L;
        FilteredData=idct(Xk); % FIBFs
        %FilteredData=real(IMFs_y);
    end
    if FigureEnable
        Xk=dct(real(FilteredData));     
        freq = Fs/2*linspace(0,1,N);
        DataFFT=abs(Xk);
        % Plot single-sided amplitude spectrum.
        %figure
        subplot(2,1,2);
        if gaindB 
        plot(freq,20*log10(DataFFT));
        else        
        plot(freq,DataFFT);
        end
        title('Filtered Data DCT plot'); 
    end
    FilteredData1=FilteredData; 
    
    
    