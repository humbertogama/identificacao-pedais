function [indmut, flag] = mutacaomodelo3(ind,percentmut)

% Funcao para mutacao de individuos em AGs.
%
%           function filmut = mutacao(fil,percentmut)
%
%           indmut     --> Individuo apos a mutacao;
%           ind        --> Individuo antes da mutacao;
%           percentmut --> Percentual de chance de ocorre mutacao em um gene;
%

tamind = length(ind);   % Determina o tamanho do indivíduo (cromossomo)
flag = 0;
if percentmut > 0
    mut = rand(1);
    if mut <= percentmut
        flag = 1;
        posmut = round(1+rand(1)*10);

        if posmut == 1
            ind(posmut) = ind(posmut) + rand(1)*2 - 1;
        elseif posmut == 2 
            ind(posmut) = abs(ind(posmut) + rand(1)*0.2 - 0.1);
        elseif posmut == 3
            ind(posmut) = ind(posmut) + rand(1)*0.2 - 0.1; 
        elseif posmut == 4
            ind(posmut) = ind(posmut) + rand(1)*0.2 - 0.1; 
        elseif posmut == 5
            ind(posmut) = ind(posmut) + rand(1)*0.2 - 0.1; 
        elseif posmut == 6
            ind(posmut) = ind(posmut) + rand(1)*2 - 1; 
        elseif posmut == 7
            ind(posmut) = ind(posmut) + rand(1)*40 - 20; 
        elseif posmut == 8
            ind(posmut) =ind(posmut) + rand(1)*600 - 300;
        elseif posmut == 9
            ind(posmut) = ind(posmut) + rand(1)*1000 - 500; 
        elseif posmut == 10 
            ind(posmut) = ind(posmut) + rand(1)*40 - 20; 
        else
            ind(posmut) = abs(ind(posmut) + rand(1)*2 - 1) ; 
        end
        
        
    end
end
indmut = ind;

end
