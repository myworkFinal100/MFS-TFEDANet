function z = sp_IDCT_FSAS(x)

    % check & convert row vector to column vector 
    [m,n] = size(x);
    if(m<n)
       x=x';
    end

    N = length(x);
    if(N == 1)
        z = x;
    else
        fx2N = fft(x, 2*N); % DCT2 using 2N FFT
        ak=(1/sqrt(N)) * [1; sqrt(2)*ones(N-1,1)]; % Normalization factor
        dctx = ak .*real( fx2N(1:N) .* exp(-1i*pi*(0:N-1)'/(2*N)) );   % DCT of x     
        
        % IDCT2 using 4N FFT
        ck=ak.*dctx;
        fx4N=fft(ck,4*N);           
        xfsas=fx4N(2:2:(2*N));% FSAS 
        
        xr=real(xfsas);   % signal x          
        xq=-imag(xfsas);  % quadrature component, - sign due to -j in exp(-j*2*pi*k*n/(N))
        z= xr+1i*xq;       
    end
    