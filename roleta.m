function sel = roleta(rol,bola)

% Funcao para selecao de pais pelo metodo da roleta em AGs.
%
%           function sel = roleta(rol,bola)
%
%           sel  = Individuo selecionado;
%           rol  = Estrutura da roleta;
%           bola = Valor sorteado;
%

flag = 0; i = 1;

while flag == 0
    if bola <= rol(i)
        sel = i; flag = 1;
    else
        i = i + 1;
    end
end
