import glob
import shutil
import os

def extract_gbff_files(source_dir, target_dir):
    #Expande o caminho com ~ para o diretório home do usuário
    source_dir = os.path.expanduser(source_dir)
    
    #Cria o diretório de destino, se ele não existir
    if not os.path.exists(target_dir):
        os.makedirs(target_dir)

    #Padrão glob para encontrar todos os arquivos .gbff em subpastas
    gbff_files = glob.glob(os.path.join(source_dir, '**', '*.gbff'), recursive=True)

    #Para cada arquivo .gbff encontrado
    for file in gbff_files:
        #Nome original do arquivo
        original_filename = os.path.basename(file)

        #Nome da pasta pai (diretório imediatamente acima do arquivo)
        parent_dir = os.path.basename(os.path.dirname(file))

        #Novo nome do arquivo: nome_da_pasta_original_nome_arquivo.gbff
        new_filename = f"{parent_dir}_{original_filename}"

        #Caminho de destino com novo nome
        target_file = os.path.join(target_dir, new_filename)
        
        try:
            #Copia o arquivo renomeado para o diretório de destino
            shutil.copy(file, target_file)
            print(f"Arquivo {original_filename} extraído e renomeado como {new_filename}")
        except Exception as e:
            print(f"Erro ao copiar o arquivo {original_filename}: {e}")

if __name__ == "__main__":
    #Caminho de origem com ~, apontando para a pasta de saída do Prokka
    source_directory = '/home/barbara/documents/RegPatternSearch/ncbi_dataset/data'  
    #Diretório onde os arquivos .gbff serão armazenados
    target_directory = './extracted_gbff'  

    #Executa a função para extrair os arquivos .gbff
    extract_gbff_files(source_directory, target_directory)
