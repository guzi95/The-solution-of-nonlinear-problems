Olá! Aqui é a Juliana Furtado, autora dos códigos aqui contidos. 

Ao adentrar nas pastas, você pode se assustar com a quantidade de arquivos; entretanto, é menos complexo do que parece. 
Atenha-se aos arquivos .f90. Ou seja, aos arquivos "program_gnewton.f90" e "module_gnewton.f90".

Utilize o compilador "gfortran" para rodar os dois arquivos .f90 que nos interessam, da seguinte forma, no terminal: 
"gfortran program_gnewton.f90 module_gnewton.f90 -o gnewton.x"

Então, você roda o executável que foi gerado:
"./gnewton.x"

-------------- OBSERVAÇÕES --------------

Existem três pastas:
 - Chutes iniciais
 - Dados obtidos
 - Lambdas

Nelas, há arquivos .dat opcionais que você pode colocar no programa para realizar testes, fazendo as devidas alterações nos códigos-fonte.

Para melhor compreender os arquivos .dat, entre no código-fonte "program_gnewton.f90", pois aí estão descritos o que cada arquivo .dat contém.


