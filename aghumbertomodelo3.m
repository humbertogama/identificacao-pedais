


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

% Determinar par�metros do AG:
%
disp(' ')
disp('Com relacao ao AG:')
disp(' ')
disp('Para inicializa��o da popula��o, temos duas op��es:')
disp(' ')
disp('     1 -> Inicializa um nova popula��o aleatoriamente;')
disp('     2 -> Continua a otimiza��o a partir de uma popula��o salva em arquivo;')
disp('        Nesse caso, precisa existir um arquivo chamado ultima_populacao.mat')
disp(' ')
initpop = -2; 
while(initpop < 1 || initpop > 2)
    initpop = input('Qual sua op��o de inicializa��o? (1 ou 2) ');
end
if initpop == 1
    %   N�mero de indiv�duos por gera��o;
    disp(' ')
    numind = input('Entre com a quantidade de indiv�duos ');
end
disp(' ')
%   N�mero m�ximo de gera��es;
maxger = input('Entre com o n�mero m�ximo de gera��es ');

%   Percentual de muta��o;
percentmut = -2; 
while(percentmut < 0 || percentmut > 10)
    disp(' ')
    pmut = input('Entre com a probabilidade percentual de um gene sofrer mutacao [0% - 10%] ');
    percentmut = pmut / 100;
end
disp(' ')
pause(.1)               % Serve para dar tempo de atualizar tela 
%=========================================================================
% Inicializar a primeira popula��o:
ger = 1;            % Primeira Gera��o
    
tampar = 10;


if initpop == 1
    disp(['Inicializando AG. Criando gera��o: ' num2str(ger) ' de ' num2str(maxger) '  . . . . . .'])
    disp(' ')
    
    for ind=1:numind
        % Cria 'numind' vetores (ordenados) para os par�metros das MFs de
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
        

        
        % cromossomo de cada indiv�duo � dado por:
        %function [out]  = compmodel(entrada,gpre, gpost,gbias,gdry,gp, gn, kp, kn, fc1, fc2)
        cromossomo(ind,:,ger) = [gpre gpost gbias gdry a b fc1 fc2 fc3 fc4 gbp];
        % =================================================================
        % =================================================================
        % || Uma solu��o mais geral para a composi��o do cromossomo      ||
        % || inicia com a defini��o autom�tica dos tamanhos dos vetores  ||
        % || de par�metros das entradas e sa�das e termina com a         ||
        % || composi��o do vetor chamado 'cromossomo'                    ||
        % =================================================================
        % =================================================================
    end
    tamcromo = length(cromossomo(1,:,1));   %  Determina o tamanho de cada 
                                            % cromossomo. Essa informa��o
                                            % ser� �til no AG.
                                            % Especialmente, no crossover.
else
    disp(['Inicializando AG. Carregando gera��o: ' num2str(ger) ' de ' num2str(maxger) '  . . . . . .'])
    disp(' ')
    load ultima_populacao;                  % Carrega uma popula��o salva
    cromossomo(:,:,ger) = ultpop;           % em arquivo
    [numind tamcromo] = size(ultpop);       % Determina o n�mero de 
                                            % indiv�duos na popula��o e o
                                            % tamanho de cada cromossomo.
    clear ultpop;                           % Apaga a vari�vel utilizada
                                            % para realizar a leitura da
                                            % popula��o salva em arquivo.
end
    
% Avaliar indiv�duos
% ------------------------------------------------------------------------
%   1. Para cada indiv�duo:


for ind = 1:numind

%       1.1. A partir do vetor Cromossomo(ind,:,ger) gerar, e escrever na
%       estrutura chamada "fuzzy_T5552", os parametros de cada MF, de cada
%       vari�vel de entrada e de sa�da).

        
% ------------------------------------------------------------------------
%       1.2. Rodar a simula��o para obter a resposta do sistema com o

    saidateste = compmodel3( x , cromossomo(ind,1,ger), cromossomo(ind,2,ger), cromossomo(ind,3,ger), cromossomo(ind,4,ger), cromossomo(ind,5,ger), cromossomo(ind,6,ger), cromossomo(ind,7,ger), cromossomo(ind,8,ger), cromossomo(ind,9,ger), cromossomo(ind,10,ger), cromossomo(ind,11,ger));

%       
% ------------------------------------------------------------------------

%       1.3. Avaliar o controlador (indiv�duo) testado.
% Obs.: Podemos usar o crit�rio de Goohart.

    %   Avalia��o de desempenho com base no Crit�rio de Goodhart para cada
    % indiv�duo, de cada gera��o. (Para detalhes, digite HELP GOODHART
    %
    % Escolha dos coeficientes (ALFAS) do �ndice de Goodharte. Esses
    % coeficientes ser�o usados, a seguir, na avalia��o dos indiv�duos da
    % 1a gera��o, e, mais adiante, na avalia��o dos indiv�duos das gera��es
    % seguintes.
    %alfa1 = 1; alfa2 = 4; alfa3 = 5; PADRAO, SINAL DE CONTROLE MUITO ALTO
    %alfa1 = 2; alfa2 = 3; alfa3 = 5; MELHOR AT� AGR
    %alfa1 = 0.2; alfa2 = 0.3; alfa3 = 1; ficou bom pra SP de 2
    
    
    
 
    alfa1 = 0.2; alfa2 = 0.3; alfa3 = 0.5;
    
    errot(ind,ger) = sqrt(sum((saidaref-saidateste).^2)/length(saidateste));
    %errot(ind,ger) = sum((saidaref-saidateste).^2)/sum(saidaref.^2);
    errofft(ind,ger) = sum((abs(fft(saidaref))-abs(fft(saidateste))).^2)/length(saidateste);
    erro(ind,ger) = pesot*errot(ind,ger) + pesofft*errofft(ind,ger);
    
    disp(['Gera��o: ' num2str(ger) ', Indiv�duo: ' num2str(ind) ' ==> Erro: ' num2str(erro(ind,ger))])
    
% ------------------------------------------------------------------------

%       1.4. Calcular a adaptabilidade o controlador (indiv�duo) testado.
% Obs.: como queremos o menor �ndice de desempenho poss�vel, a adap. de
% cada indiv�duo ser� o inverso da fun��o de avalia��o:

    adap(ind,ger) = 1 / erro(ind,ger);    % Calculo do �ndice de Goodhart
                                        % para cada indiv�duo, de cada
                                        %gera��o.

end
% ------------------------------------------------------------------------
%   2. Calcular a adaptabilidade relativa dos indiv�duos dessa gera��o:

adaprel(:,ger) = adap(:,ger)/sum(adap(:,ger));

% ------------------------------------------------------------------------
%   3. Determina quem s�o o melhor e o pior indiv�duo dessa gera��o:

[melhoradap(ger),posmelhor(ger)] = max(adap(:,ger));    % Determina o valor
                                                        % da melhor adap.
                                                        % da gera��o, e a
                                                        % posi��o do ind.
                                                        % com tal adap.
                                                        % dentro da gera��o
melhorind(ger,:) = cromossomo(posmelhor(ger),:,ger);    % Garda o individuo
                                                        % mais adaptado de
                                                        % cada gera��o

[pioradap(ger),pospior(ger)] = min(adap(:,ger));        % Determina o valor
                                                        % da pior adap.
                                                        % da gera��o, e a
                                                        % posi��o do ind.
                                                        % com tal adap.
                                                        % dentro da gera��o
piorind(ger,:) = cromossomo(pospior(ger),:,ger);        % Garda o individuo
                                                        % menos adaptado de
                                                        % cada gera��o

% Fim da primeira popula��o
%=========================================================================

mediapop(ger) = sum(adap(:,ger))/numind;	% Calcula uma m�dia das
                                            % adaptabilidades da gera��o
                                            %
                                            % Obs. Podemos pensar em usar
                                            % um arredondamento em 'n'
                                            % casas decimais.
                                           
% Cria vari�veis para alguns dos crit�rios de parada
mediapopigual = 0; mesmo = 0;


%=========================================================================
%||                     In�cio do Loop do AG                            ||
%=========================================================================
% Gerar popula��es seguintes
while (mediapopigual <= 400 & mesmo <= 100 & ger < maxger)   % crit�rios de
    % parada para o AG:
    %   - Estagna��o da m�dia da popula��o por, pelo menos, 40 gera��es;
    %   - Estagna��o do melhor indiv�duo por, pelo menos, 80 gera��es;
    %   - Atingir o n�mero m�ximo de gera��es determinado pelo usu�rio;

    disp(' ')
    disp(['Rodando AG. Criando gera��o: ' num2str(ger+1) ' de ' num2str(maxger) '  . . . . . . . .'])
    disp(' ')
    
    for i = 1:numind                        % Cria a estrtura da roleta
        rol(i,ger) = sum(adaprel(1:i,ger)); % Soma a adap. rel. de cada
                                            % indiv�duo coma a adap.rel.
                                            % dos anteriores, para
                                            % determinar a por��o da roleta
                                            % que cabe a cada indiv�duo. O
                                            % primeiro ind. vai de zero at�
                                            % sua adap.rel.. O segundo, vai
                                            % da adap.rel. do primeiro at�
                                            % a soma da adap.rel do
                                            % primeiro com a sua. E assim
                                            % por diante.
    end
    
   
    for ind = 1:2:2*floor(numind/2)        % faz 'ind' variar de 2 em 2, do 
                                           % primeiro indiv�duo at� o
                                           % pen�timo indiv�duo, se o total
                                           % de indiv�duos for par, ou at�
                                           % o anti-pen�ltimo, se o total
                                           % for impar
                                           % Se for impar, vai ficar
                                           % faltando um indiv�duo para
                                           % completar a nova popula��o.
                                           % Uma possibilidade � clonar o
                                           % melhor indiv�duo da popula��o
                                           % atual.
                                           
       bola(ind:ind+1,ger) = rand(1,2);    % 'Arremessa duas bolas sobre a
                                           % roleta, para sortear dois
                                           % indiv�duos que ir�o gerar dois
                                           % descendentes para pr�xima
                                           % gera��o.
                                           
    
       poscros(ind,ger) = round(10*rand(1) + 1); %1 a 11
       % Gera o primeiro descendente do crossover
       cromossomo(ind,:,ger+1) = [cromossomo(roleta(rol(:,ger),bola(ind,ger)),1:poscros(ind,ger),ger) cromossomo(roleta(rol(:,ger),bola(ind+1,ger)),poscros(ind,ger)+1:tamcromo,ger)];
       % Mutacao no primeiro descendente do crossover
       [cromossomo(ind,:,ger+1), mflag] = mutacaomodelo3(cromossomo(ind,:,ger+1),percentmut);
       if mflag == 1
           disp(['Houve mutacao no individuo: ', num2str(ind),' da gera��o (descendencia): ', num2str(ger+1)]);
       end
       % Gera o segundo descendente do crossover
       cromossomo(ind+1,:,ger+1) = [cromossomo(roleta(rol(:,ger),bola(ind+1,ger)),1:poscros(ind,ger),ger) cromossomo(roleta(rol(:,ger),bola(ind,ger)),poscros(ind,ger)+1:tamcromo,ger)];
       % Mutacao no segundo descendente do crossover
       [cromossomo(ind+1,:,ger+1), mflag1] = mutacaomodelo3(cromossomo(ind+1,:,ger+1),percentmut);
       if mflag1 == 1
           disp(['Houve mutacao no individuo: ', num2str(ind+1),' da gera��o (descendencia): ', num2str(ger+1)]);
       end    
       
    end
    
    % Clonagem do melhor indiv�duo
    if round(numind/2) == numind/2                      % Se o n�mero total
                                                        % de indiv�duos da
                                                        % popula��o for par
       indclone = round(rand(1) * (numind - 1) + 1);    % Substitui um ind.
                                                        % aleatorio da nova
                                                        % popula��o, pelo
       cromossomo(indclone,:,ger+1) = melhorind(ger,:); % melho ind. da
                                                        % pop. anterior.
    else                                                % Se o n�mero total
                                                        % de indiv�duos da
                                                        % pop. for impar
       cromossomo(numind,:,ger+1) = melhorind(ger,:);   % Completa a pop. 
                                                        % com o melhor ind.
                                                        % da pop. anterior
    end
    
    
    
    % Avaliar indiv�duos
    % ---------------------------------------------------------------------
    %   1. Para cada indiv�duo da nova gera��o (ger + 1):

    for ind = 1:numind
        


%       1.1. A partir do vetor Cromossomo(ind,:,ger+1) gerar, e escrever na
%       estrutura chamada "fuzzy_T5552", os parametros de cada MF, de cada
%       vari�vel de entrada e de sa�da).

        % -----------------------------------------------------------------
        %       1.2. Rodar a simula��o para obter a resposta do sistema com
        %  o controlador (indiv�duo) a ser avaliado.

        saidateste = compmodel3( x , cromossomo(ind,1,ger+1), cromossomo(ind,2,ger+1), cromossomo(ind,3,ger+1), cromossomo(ind,4,ger+1), cromossomo(ind,5,ger+1), cromossomo(ind,6,ger+1), cromossomo(ind,7,ger+1), cromossomo(ind,8,ger+1), cromossomo(ind,9,ger+1), cromossomo(ind,10,ger+1), cromossomo(ind,11,ger+1));

    
        % -----------------------------------------------------------------
        %       1.3. Avaliar o controlador (indiv�duo) testado.
        % Avalia��o de desempenho com base no Crit�rio de Goodhart para
        % cada indiv�duo da nova gera��o(ger+1).

       errot(ind,ger+1) = sqrt(sum((saidaref-saidateste).^2)/length(saidateste));
        %errot(ind,ger+1) = sum((saidaref-saidateste).^2)/sum(saidaref.^2);
        errofft(ind,ger+1) = sum((abs(fft(saidaref))-abs(fft(saidateste))).^2)/length(saidateste);
        erro(ind,ger+1) = pesot*errot(ind,ger+1) + pesofft*errofft(ind,ger+1);

        disp(['Gera��o: ' num2str(ger+1) ', Indiv�duo: ' num2str(ind) ' ==> Erro: ' num2str(erro(ind,ger+1))])
        
        % -----------------------------------------------------------------
        %       1.4. Calcular a adaptabilidade do controlador (indiv�duo).
    
        adap(ind,ger+1) = 1 / erro(ind,ger+1);    % Calculo da adaptabilidade
                                                % de cada indiv�duo
    end
    % ---------------------------------------------------------------------
    %   2. Calcular a adaptabilidade relativa dos indiv�duos dessa gera��o:

    adaprel(:,ger+1) = adap(:,ger+1)/sum(adap(:,ger+1));% Calculo da 
                                                        % adaptabilidade 
                                                        % relativa de cada
                                                        % indiv�duo

    % ---------------------------------------------------------------------
    %   3. Determina quem � o melhor indiv�duo dessa gera��o:

    [melhoradap(ger+1),posmelhor(ger+1)] = max(adap(:,ger+1));  % Determina
    % o valor da melhor adap. da gera��o, e a posi��o do ind. com tal adap.
    % dentro da gera��o
    melhorind(ger+1,:) = cromossomo(posmelhor(ger+1),:,ger+1);	% Garda o 
    % individuo mais adaptado da gera��o

    % Fim da n-esima popula��o
    %======================================================================

    % Calcula vari�veis para alguns dos crit�rios de parada
    mediapop(ger+1) = sum(adap(:,ger+1))/numind;% Calcula uma m�dia das
                                                % adaptabilidades da ger.
    % Avalia a estagna��o da Popula��o
	if mediapop(ger+1) == mediapop(ger)
        mediapopigual = mediapopigual + 1;      
    else
        mediapopigual = 0;
    end
    % Avalia a estagna��o do melhor indiv�duo
    if melhorind(ger) == melhorind(ger+1)
        mesmo = mesmo + 1;
    else
        mesmo = 0;
    end 
    ger = ger + 1;  % Incrementa o indexador de gera��o, para, se nenhum
                    % crit�rio de parada for atingido, retornar ao in�cio
                    % do loop (while) e criar uma nova popula��o.              
    %======================================================================
end
%=========================================================================
%||                        Fim do Loop do AG                            ||
%=========================================================================
disp('Melhor indiv�duo: ')
melhorind(ger,:)
                                             
% Salva a �ltima popula��o em arquivo, para utiliza��es posteriores, com o
% nome de 'ultima_populacao'.

ultpop = cromossomo(:,:,ger);           % �ltima popula��o gerada pelo AG
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
disp('Par�metros: ')
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

