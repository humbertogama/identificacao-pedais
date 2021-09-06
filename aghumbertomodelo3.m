


clear all
close all
clc




pesot = 0;
pesofft = 1;

% [x, fs] = audioread('cleanr.wav');
% [saidaref, fs] = audioread('jmpr.wav');

[x, fs] = audioread('testeBDL.wav');
[saidaref, fs] = audioread('testeBDS.wav');
entrada = x;

% Determinar parâmetros do AG:
%
disp(' ')
disp('Com relacao ao AG:')
disp(' ')
disp('Para inicialização da população, temos duas opções:')
disp(' ')
disp('     1 -> Inicializa um nova população aleatoriamente;')
disp('     2 -> Continua a otimização a partir de uma população salva em arquivo;')
disp('        Nesse caso, precisa existir um arquivo chamado ultima_populacao.mat')
disp(' ')
initpop = -2; 
while(initpop < 1 || initpop > 2)
    initpop = input('Qual sua opção de inicialização? (1 ou 2) ');
end
if initpop == 1
    %   Número de indivíduos por geração;
    disp(' ')
    numind = input('Entre com a quantidade de indivíduos ');
end
disp(' ')
%   Número máximo de gerações;
maxger = input('Entre com o número máximo de gerações ');

%   Percentual de mutação;
percentmut = -2; 
while(percentmut < 0 || percentmut > 10)
    disp(' ')
    pmut = input('Entre com a probabilidade percentual de um gene sofrer mutacao [0% - 10%] ');
    percentmut = pmut / 100;
end
disp(' ')
pause(.1)               % Serve para dar tempo de atualizar tela 
%=========================================================================
% Inicializar a primeira população:
ger = 1;            % Primeira Geração
    
tampar = 10;


if initpop == 1
    disp(['Inicializando AG. Criando geração: ' num2str(ger) ' de ' num2str(maxger) '  . . . . . .'])
    disp(' ')
    
    for ind=1:numind
        % Cria 'numind' vetores (ordenados) para os parâmetros das MFs de
        % cada entrada
        
        %compmodel3(entrada,gpre, gpost,gbias,gdry,a, b, fc1, fc2, fc3, fc4, gbp)
        
      
        gpre = 1+rand(1)*9; %entre 1 e 10
        gpost = 0.2+rand(1)*0.8; %entre 0.2 e 1
        gbias = rand(1); %entre 0 e 1
        gdry = rand(1); %entre 0 e 1
        a = 0.2+0.8*rand(1); %entre 0.2 e 1
        b = 1+rand(1)*10; %entre 1 e 11
        fc1 = 60 + rand(1)*140; %entre 60 e 200 Hz
        fc2 = 500+rand(1)*2500; %entre 500 e 3000 Hz
        fc3 = 1000 + rand(1)*4000; %entre 1000 e 5000 Hz
        fc4 = 30+rand(1)*170; %entre 30 e 200 Hz
        gbp = rand(1)*15; %entre 1 e 15;
        

        
        % cromossomo de cada indivíduo é dado por:
        %function [out]  = compmodel(entrada,gpre, gpost,gbias,gdry,gp, gn, kp, kn, fc1, fc2)
        cromossomo(ind,:,ger) = [gpre gpost gbias gdry a b fc1 fc2 fc3 fc4 gbp];
        % =================================================================
        % =================================================================
        % || Uma solução mais geral para a composição do cromossomo      ||
        % || inicia com a definição automática dos tamanhos dos vetores  ||
        % || de parâmetros das entradas e saídas e termina com a         ||
        % || composição do vetor chamado 'cromossomo'                    ||
        % =================================================================
        % =================================================================
    end
    tamcromo = length(cromossomo(1,:,1));   %  Determina o tamanho de cada 
                                            % cromossomo. Essa informação
                                            % será útil no AG.
                                            % Especialmente, no crossover.
else
    disp(['Inicializando AG. Carregando geração: ' num2str(ger) ' de ' num2str(maxger) '  . . . . . .'])
    disp(' ')
    load ultima_populacao;                  % Carrega uma população salva
    cromossomo(:,:,ger) = ultpop;           % em arquivo
    [numind tamcromo] = size(ultpop);       % Determina o número de 
                                            % indivíduos na população e o
                                            % tamanho de cada cromossomo.
    clear ultpop;                           % Apaga a variável utilizada
                                            % para realizar a leitura da
                                            % população salva em arquivo.
end
    
% Avaliar indivíduos
% ------------------------------------------------------------------------
%   1. Para cada indivíduo:


for ind = 1:numind

%       1.1. A partir do vetor Cromossomo(ind,:,ger) gerar, e escrever na
%       estrutura chamada "fuzzy_T5552", os parametros de cada MF, de cada
%       variável de entrada e de saída).

        
% ------------------------------------------------------------------------
%       1.2. Rodar a simulação para obter a resposta do sistema com o

    saidateste = compmodel3( x , cromossomo(ind,1,ger), cromossomo(ind,2,ger), cromossomo(ind,3,ger), cromossomo(ind,4,ger), cromossomo(ind,5,ger), cromossomo(ind,6,ger), cromossomo(ind,7,ger), cromossomo(ind,8,ger), cromossomo(ind,9,ger), cromossomo(ind,10,ger), cromossomo(ind,11,ger));

%       
% ------------------------------------------------------------------------

%       1.3. Avaliar o controlador (indivíduo) testado.
% Obs.: Podemos usar o critério de Goohart.

    %   Avaliação de desempenho com base no Critério de Goodhart para cada
    % indivíduo, de cada geração. (Para detalhes, digite HELP GOODHART
    %
    % Escolha dos coeficientes (ALFAS) do índice de Goodharte. Esses
    % coeficientes serão usados, a seguir, na avaliação dos indivíduos da
    % 1a geração, e, mais adiante, na avaliação dos indivíduos das gerações
    % seguintes.
    %alfa1 = 1; alfa2 = 4; alfa3 = 5; PADRAO, SINAL DE CONTROLE MUITO ALTO
    %alfa1 = 2; alfa2 = 3; alfa3 = 5; MELHOR ATÉ AGR
    %alfa1 = 0.2; alfa2 = 0.3; alfa3 = 1; ficou bom pra SP de 2
    
    
    
 
    alfa1 = 0.2; alfa2 = 0.3; alfa3 = 0.5;
    
    errot(ind,ger) = sqrt(sum((saidaref-saidateste).^2)/length(saidateste));
    %errot(ind,ger) = sum((saidaref-saidateste).^2)/sum(saidaref.^2);
    errofft(ind,ger) = sum((abs(fft(saidaref))-abs(fft(saidateste))).^2)/length(saidateste);
    erro(ind,ger) = pesot*errot(ind,ger) + pesofft*errofft(ind,ger);
    
    disp(['Geração: ' num2str(ger) ', Indivíduo: ' num2str(ind) ' ==> Erro: ' num2str(erro(ind,ger))])
    
% ------------------------------------------------------------------------

%       1.4. Calcular a adaptabilidade o controlador (indivíduo) testado.
% Obs.: como queremos o menor índice de desempenho possível, a adap. de
% cada indivíduo será o inverso da função de avaliação:

    adap(ind,ger) = 1 / erro(ind,ger);    % Calculo do ìndice de Goodhart
                                        % para cada indivíduo, de cada
                                        %geração.

end
% ------------------------------------------------------------------------
%   2. Calcular a adaptabilidade relativa dos indivíduos dessa geração:

adaprel(:,ger) = adap(:,ger)/sum(adap(:,ger));

% ------------------------------------------------------------------------
%   3. Determina quem são o melhor e o pior indivíduo dessa geração:

[melhoradap(ger),posmelhor(ger)] = max(adap(:,ger));    % Determina o valor
                                                        % da melhor adap.
                                                        % da geração, e a
                                                        % posição do ind.
                                                        % com tal adap.
                                                        % dentro da geração
melhorind(ger,:) = cromossomo(posmelhor(ger),:,ger);    % Garda o individuo
                                                        % mais adaptado de
                                                        % cada geração

[pioradap(ger),pospior(ger)] = min(adap(:,ger));        % Determina o valor
                                                        % da pior adap.
                                                        % da geração, e a
                                                        % posição do ind.
                                                        % com tal adap.
                                                        % dentro da geração
piorind(ger,:) = cromossomo(pospior(ger),:,ger);        % Garda o individuo
                                                        % menos adaptado de
                                                        % cada geração

% Fim da primeira população
%=========================================================================

mediapop(ger) = sum(adap(:,ger))/numind;	% Calcula uma média das
                                            % adaptabilidades da geração
                                            %
                                            % Obs. Podemos pensar em usar
                                            % um arredondamento em 'n'
                                            % casas decimais.
                                           
% Cria variáveis para alguns dos critérios de parada
mediapopigual = 0; mesmo = 0;


%=========================================================================
%||                     Início do Loop do AG                            ||
%=========================================================================
% Gerar populações seguintes
while (mediapopigual <= 400 & mesmo <= 100 & ger < maxger)   % critérios de
    % parada para o AG:
    %   - Estagnação da média da população por, pelo menos, 40 gerações;
    %   - Estagnação do melhor indivíduo por, pelo menos, 80 gerações;
    %   - Atingir o número máximo de gerações determinado pelo usuário;

    disp(' ')
    disp(['Rodando AG. Criando geração: ' num2str(ger+1) ' de ' num2str(maxger) '  . . . . . . . .'])
    disp(' ')
    
    for i = 1:numind                        % Cria a estrtura da roleta
        rol(i,ger) = sum(adaprel(1:i,ger)); % Soma a adap. rel. de cada
                                            % indivíduo coma a adap.rel.
                                            % dos anteriores, para
                                            % determinar a porção da roleta
                                            % que cabe a cada indivíduo. O
                                            % primeiro ind. vai de zero até
                                            % sua adap.rel.. O segundo, vai
                                            % da adap.rel. do primeiro até
                                            % a soma da adap.rel do
                                            % primeiro com a sua. E assim
                                            % por diante.
    end
    
   
    for ind = 1:2:2*floor(numind/2)        % faz 'ind' variar de 2 em 2, do 
                                           % primeiro indivíduo até o
                                           % penútimo indivíduo, se o total
                                           % de indivíduos for par, ou até
                                           % o anti-penúltimo, se o total
                                           % for impar
                                           % Se for impar, vai ficar
                                           % faltando um indivíduo para
                                           % completar a nova população.
                                           % Uma possibilidade é clonar o
                                           % melhor indivíduo da população
                                           % atual.
                                           
       bola(ind:ind+1,ger) = rand(1,2);    % 'Arremessa duas bolas sobre a
                                           % roleta, para sortear dois
                                           % indivíduos que irão gerar dois
                                           % descendentes para próxima
                                           % geração.
                                           
    
       poscros(ind,ger) = round(10*rand(1) + 1); %1 a 11
       % Gera o primeiro descendente do crossover
       cromossomo(ind,:,ger+1) = [cromossomo(roleta(rol(:,ger),bola(ind,ger)),1:poscros(ind,ger),ger) cromossomo(roleta(rol(:,ger),bola(ind+1,ger)),poscros(ind,ger)+1:tamcromo,ger)];
       % Mutacao no primeiro descendente do crossover
       [cromossomo(ind,:,ger+1), mflag] = mutacaomodelo3(cromossomo(ind,:,ger+1),percentmut);
       if mflag == 1
           disp(['Houve mutacao no individuo: ', num2str(ind),' da geração (descendencia): ', num2str(ger+1)]);
       end
       % Gera o segundo descendente do crossover
       cromossomo(ind+1,:,ger+1) = [cromossomo(roleta(rol(:,ger),bola(ind+1,ger)),1:poscros(ind,ger),ger) cromossomo(roleta(rol(:,ger),bola(ind,ger)),poscros(ind,ger)+1:tamcromo,ger)];
       % Mutacao no segundo descendente do crossover
       [cromossomo(ind+1,:,ger+1), mflag1] = mutacaomodelo3(cromossomo(ind+1,:,ger+1),percentmut);
       if mflag1 == 1
           disp(['Houve mutacao no individuo: ', num2str(ind+1),' da geração (descendencia): ', num2str(ger+1)]);
       end    
       
    end
    
    % Clonagem do melhor indivíduo
    if round(numind/2) == numind/2                      % Se o número total
                                                        % de indivíduos da
                                                        % população for par
       indclone = round(rand(1) * (numind - 1) + 1);    % Substitui um ind.
                                                        % aleatorio da nova
                                                        % população, pelo
       cromossomo(indclone,:,ger+1) = melhorind(ger,:); % melho ind. da
                                                        % pop. anterior.
    else                                                % Se o número total
                                                        % de indivíduos da
                                                        % pop. for impar
       cromossomo(numind,:,ger+1) = melhorind(ger,:);   % Completa a pop. 
                                                        % com o melhor ind.
                                                        % da pop. anterior
    end
    
    
    
    % Avaliar indivíduos
    % ---------------------------------------------------------------------
    %   1. Para cada indivíduo da nova geração (ger + 1):

    for ind = 1:numind
        


%       1.1. A partir do vetor Cromossomo(ind,:,ger+1) gerar, e escrever na
%       estrutura chamada "fuzzy_T5552", os parametros de cada MF, de cada
%       variável de entrada e de saída).

        % -----------------------------------------------------------------
        %       1.2. Rodar a simulação para obter a resposta do sistema com
        %  o controlador (indivíduo) a ser avaliado.

        saidateste = compmodel3( x , cromossomo(ind,1,ger+1), cromossomo(ind,2,ger+1), cromossomo(ind,3,ger+1), cromossomo(ind,4,ger+1), cromossomo(ind,5,ger+1), cromossomo(ind,6,ger+1), cromossomo(ind,7,ger+1), cromossomo(ind,8,ger+1), cromossomo(ind,9,ger+1), cromossomo(ind,10,ger+1), cromossomo(ind,11,ger+1));

    
        % -----------------------------------------------------------------
        %       1.3. Avaliar o controlador (indivíduo) testado.
        % Avaliação de desempenho com base no Critério de Goodhart para
        % cada indivíduo da nova geração(ger+1).

       errot(ind,ger+1) = sqrt(sum((saidaref-saidateste).^2)/length(saidateste));
        %errot(ind,ger+1) = sum((saidaref-saidateste).^2)/sum(saidaref.^2);
        errofft(ind,ger+1) = sum((abs(fft(saidaref))-abs(fft(saidateste))).^2)/length(saidateste);
        erro(ind,ger+1) = pesot*errot(ind,ger+1) + pesofft*errofft(ind,ger+1);

        disp(['Geração: ' num2str(ger+1) ', Indivíduo: ' num2str(ind) ' ==> Erro: ' num2str(erro(ind,ger+1))])
        
        % -----------------------------------------------------------------
        %       1.4. Calcular a adaptabilidade do controlador (indivíduo).
    
        adap(ind,ger+1) = 1 / erro(ind,ger+1);    % Calculo da adaptabilidade
                                                % de cada indivíduo
    end
    % ---------------------------------------------------------------------
    %   2. Calcular a adaptabilidade relativa dos indivíduos dessa geração:

    adaprel(:,ger+1) = adap(:,ger+1)/sum(adap(:,ger+1));% Calculo da 
                                                        % adaptabilidade 
                                                        % relativa de cada
                                                        % indivíduo

    % ---------------------------------------------------------------------
    %   3. Determina quem é o melhor indivíduo dessa geração:

    [melhoradap(ger+1),posmelhor(ger+1)] = max(adap(:,ger+1));  % Determina
    % o valor da melhor adap. da geração, e a posição do ind. com tal adap.
    % dentro da geração
    melhorind(ger+1,:) = cromossomo(posmelhor(ger+1),:,ger+1);	% Garda o 
    % individuo mais adaptado da geração

    % Fim da n-esima população
    %======================================================================

    % Calcula variáveis para alguns dos critérios de parada
    mediapop(ger+1) = sum(adap(:,ger+1))/numind;% Calcula uma média das
                                                % adaptabilidades da ger.
    % Avalia a estagnação da População
	if mediapop(ger+1) == mediapop(ger)
        mediapopigual = mediapopigual + 1;      
    else
        mediapopigual = 0;
    end
    % Avalia a estagnação do melhor indivíduo
    if melhorind(ger) == melhorind(ger+1)
        mesmo = mesmo + 1;
    else
        mesmo = 0;
    end 
    ger = ger + 1;  % Incrementa o indexador de geração, para, se nenhum
                    % critério de parada for atingido, retornar ao início
                    % do loop (while) e criar uma nova população.              
    %======================================================================
end
%=========================================================================
%||                        Fim do Loop do AG                            ||
%=========================================================================
disp('Melhor indivíduo: ')
melhorind(ger,:)
                                             
% Salva a última população em arquivo, para utilizações posteriores, com o
% nome de 'ultima_populacao'.

ultpop = cromossomo(:,:,ger);           % Última população gerada pelo AG
save ultima_populacao ultpop


t = -1:0.001:1;
out  = saturacaocompleta( t , melhorind(ger,1), melhorind(ger,2), melhorind(ger,3), melhorind(ger,4), melhorind(ger,5), melhorind(ger,6));
out2 = saturacaocompleta( x , melhorind(ger,1), melhorind(ger,2), melhorind(ger,3), melhorind(ger,4), melhorind(ger,5), melhorind(ger,6));
outs = compmodel3(x, melhorind(ger,1), melhorind(ger,2), melhorind(ger,3), melhorind(ger,4), melhorind(ger,5), melhorind(ger,6), melhorind(ger,7), melhorind(ger,8), melhorind(ger,9), melhorind(ger,10), melhorind(ger,11));
figure(3), plot(t,out);
figure(4), plot(x,out2)

% sound(outs, fs);
% audiowrite('nome.wav', outs, 44100);


disp('RESUMO: ');
disp('Parâmetros: ')
disp(melhorind(ger,1))
disp(melhorind(ger,2))
disp(melhorind(ger,3))
disp(melhorind(ger,4))
disp(melhorind(ger,5))
disp(melhorind(ger,6))
disp(melhorind(ger,7))
disp(melhorind(ger,8))
disp(melhorind(ger,9))
disp(melhorind(ger,10))
disp(melhorind(ger,11))



% sound(outs, fs);
% audiowrite('nome.wav', outs, 44100);

