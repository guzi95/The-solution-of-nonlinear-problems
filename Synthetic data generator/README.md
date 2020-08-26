Olá! Aqui é a Juliana Furtado, autora dos códigos aqui contidos. 

Aqui os dados gravimétricos sintéticos vão ser gerados. 
Trata-se de um problema direto para modelar uma anomalia geológica esférica na subsuperfície terrestre.
Os dados que são gerados aqui, foram utilizados para testar os métodos de inversão não-linear 
Levenberg-Marquardt, Máximo Declive e Gauss-Newton.

Na pasta "Parametros para modelagem" há arquivos .dat com diferentes valores de raio e profundidade
para gerar diferentes anomalias gravimétricas. Basta trazer estes arquivos para o diretório
principal "Synthetic data generator" e inserir o nome desses arquivos na linha de código:

 - open (file='parametros_a.dat', status='old', unit=10, action='read')
 
 , no lugar de "parametros_a.dat"

Na pasta "Dados gerados" estão os valores de anomalias gravimétricas geradas por estes testes.

