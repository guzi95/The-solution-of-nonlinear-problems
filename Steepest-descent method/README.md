Olá! Aqui é a Juliana Furtado, autora dos códigos aqui contidos. 

Ao adentrar nas pastas, você pode se assustar com a quantidade de arquivos; 
entretanto, é menos complexo do que parece. 
Atenha-se aos arquivos .f90. Ou seja, aos arquivos "program_steepest_descf90" e "module_steepest_desc.f90".

Utilize o compilador "gfortran" para rodar no terminal os dois arquivos .f90 que nos interessam, da seguinte forma: 

 - gfortran program_steepest_desc.f90 module_steepest_desc.f90 -o steepest.x

Então, você roda o executável que foi gerado:
 - ./steepest.x
 
 O resultado final encontrado é imprimido no final do arquivo "final_x.dat".
 Nele, estão contidos as estimativas obtidas pelo método do Máximo Declive durante toda a iteração.

-------------- OBSERVAÇÕES --------------

Existem três pastas:
 - Chutes iniciais
 - Dados obtidos
 - Lambdas
 - Valores de u (leia "mi")

Nelas, há arquivos .dat opcionais que você pode colocar no programa para realizar testes, fazendo as devidas alterações nos códigos-fonte.

Para melhor compreender os arquivos .dat, entre no código-fonte "program_steepest_desc.f90", pois aí estão descritos o que cada arquivo .dat contém.




