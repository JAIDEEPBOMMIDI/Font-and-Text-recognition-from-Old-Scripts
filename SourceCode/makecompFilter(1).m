function filt=makecompFilter(N)
% Function to generate coefficients of Daubechies symmetric Complex orthonormal wavelet filters.  
% Input: N -- should be of the form 2 + 4*m, % m=1,2,3,...
% e-mail: issacniwas@gmail.com
    if mod(N,2) ~=0
       error('N must be even');
    end
    
    p=N/2;  %works up to p=80
    b(p)=1;
    for i=p-1:-1:1
       b(i)= b(i+1)*(2*p-i-1)/(4*(p-i));
    end

    r=roots(b)/4;  
    %convert from y to z:
    z1= (1-2*r) + sqrt((1-2*r).^2 -1);
    z2= (1-2*r) - sqrt((1-2*r).^2 -1);  
    z=[z1;z2];
    zchoose = z( abs(z)>1 & imag(z) >0);
    [junk,index]= sort(real(zchoose));   
    zchoose = zchoose(index);           
    zchoose(1:2:length(zchoose))=conj(zchoose(1:2:length(zchoose)));
    zchoose = [zchoose; 1./zchoose];
    if(length(zchoose) ~= length(z)/2)
       error('N not of form 2+4m : No symmetric solutions for this N')
    end
    q=poly(zchoose);
    %add p zeros at -1
    flat=1;
    for i=1:p
       flat=conv(flat,[1 1]);
    end
    dau= conv(q,flat);
    filt= sqrt(2)*dau/sum(dau);