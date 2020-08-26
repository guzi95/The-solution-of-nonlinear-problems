Olá! Aqui é a Juliana Furtado, autora dos códigos aqui contidos. 

Ao adentrar nas pastas, você pode se assustar com a quantidade de arquivos; 
entretanto, é menos complexo do que parece. 
Atenha-se aos arquivos .f90. Ou seja, aos arquivos "program_levenberg_marquardt.f90" e "module_levenberg_marquardt.f90".

Existem três pastas:
 - Chutes iniciais
 - Dados obtidos
 - Lambdas

Nelas, há arquivos .dat opcionais que você pode colocar no programa para realizar testes, fazendo as devidas alterações nos códigos-fonte, alterando o nome dos arquivos.

Para melhor compreender os arquivos .dat nas pastas dos programas, entre no código-fonte "program_levenberg_marquardt.f90", pois aí estão descritos o que cada arquivo .dat salva.

Utilize o compilador "gfortran" para rodar os dois arquivos .f90 que nos interessam, da seguinte forma, no terminal: 

 - "gfortran program_levenberg_marquardt.f90 module_levenberg_marquardt.f90 -o lev_marq.x"

Então, você roda o executável que foi gerado:
 - "./lev_marq.x"
 
 O resultado final encontrado é imprimido no arquivo "final_x.dat".
 Neste arquivo, estão contidos as estimativas obtidas pelo método de Levenberg-Marquardt.


