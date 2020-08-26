Olá! Aqui é a Juliana Furtado, autora dos códigos aqui contidos. 

Ao adentrar nas pastas, você pode se assustar com a quantidade de arquivos; entretanto, é menos complexo do que parece. 
Atenha-se aos arquivos .f90. Ou seja, as arquivos "program_steepest_desc.f90" e "module_steepest_desc.f90".

Existem quatro pastas:
 - Chutes iniciais
 - Dados obtidos
 - Lambdas
 - Valores de mi

Nelas, há arquivos .dat opcionais que você pode rodar no programa para realizar testes, fazendo as devidas alterações nos códigos-fonte, alterando o nome dos arquivos de entrada e saída que você quer utilizar.

Para compreender os arquivos .dat nas pastas dos programas, entre no código-fonte "program_steepest_desc.f90", pois aí estão descritos o que cada arquivo .dat salva ou abre.

Utilize o compilador "gfortran" para rodar os dois arquivos .f90 que nos interessam, da seguinte forma, no terminal: 
"gfortran program_steepest_desc.f90 module_steepest_desc.f90 -o steepest.x"

Então, você roda o executável que foi gerado:
"./steepest.x"


