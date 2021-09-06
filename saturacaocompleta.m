function [out]  = saturacaocompleta(x,gpre, gpost,gbias,gdry,a,b)
    gbias = 0;
    gdry = 0;
    fs = 44100;
    x = x*gpre;
    gwet = 1-gdry;
    
    
    lowpass = fir1(50,5000/(fs/2));
    xabs = abs(x);
    xlp = filter(lowpass,1,xabs);
    xlp = gbias*xlp;
    xmap = x-xlp;
    
    %mapeamento
    mapout = xmap;
    tam = length(x);
    for i = 1:tam
        mapout(i) = a*tanh(b*xmap(i)); 
    end   
    
    mapout = gwet*mapout;
    out = mapout + gdry*x;
    out = out*gpost;
end

