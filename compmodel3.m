function [out]  = compmodel3(entrada,gpre, gpost,gbias,gdry,a,b, fc1, fc2, fc3, fc4, gbp)
    fs = 44100;
    np = 6;
    %estrutura: BPF -> NL -> LPF -> HPF
    bandpass1 = fir1(np,[fc1/22050 fc2/22050],'DC-0');
    lowpass2 = fir1(np,fc3/(fs/2));
    highpass2 = fir1(np, fc4/(fs/2),'high');

    gdry = 0;
    gbias = 0;

    xm = filter(bandpass1,1,entrada);
    x = entrada + gbp*xm;
    
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
     
%      %mapeamento
%     x = xmap;
%     y = x;
%     
%     for i = 1:length(x)
%         if abs(x(i)) > 1
%             x(i) = sign(x(i));
%         end    
%     end    
%     
%     for i = 1:length(x)
%         if (x(i) >= -1) && (x(i) < -0.08905)
%             y(i) = -0.75*(1-(1-(abs(x(i))-0.032847))^(12) + (1/3)*(abs(x(i))-0.032847))+0.01;
%         elseif (x(i) >= -0.08905) && (x(i) < 0.320018)
%             y(i) = -6.153*x(i)^2 + 3.9375*x(i);
%         else
%             y(i) = 0.630065;  
%         end    
%     
%     end 
%     
%     
%     mapout = gwet*y;
   
    out1 = mapout + gdry*x;
    out1 = out1*gpost;
    
    out2 = filter(lowpass2, 1, out1);
    out = filter(highpass2, 1, out2);
    
    
end

