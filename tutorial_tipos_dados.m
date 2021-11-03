%função para visualizar todas as variaveis criadas
whos

%strings sao consideradas vetores de linha única

%structs sao similares a dicionarios
my_struct.name = 'My new struct'
class(my_struct)

mystruct.age = 25

%diz se o label é uma informação da struct 
isfield(my_struct, 'name')
%remove uma label do struct
rmfield(my_struct, 'age')

%para poder definir um campo on the go
setfield(my_struct,'gemder','f')

%struct dentro de struct
my_struct.contact.phone = 12345567890
my_struct.contact.email = 'person@gmail.com'

%declaracao multipla de campos de uma struct

S = struct('nome','luan','email','luan@gmail.com')

%cells uteis para armazenar coisas pelo nomee coisas muito distintas
my_cell{1} = 'hello world!';
my_cell{'A'} = [1 2; 3 4];

