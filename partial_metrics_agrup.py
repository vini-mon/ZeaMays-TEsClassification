import os
import pandas as pd
from tqdm import tqdm

def agrupar_csvs(pasta, save_path):
    
    if not os.path.exists(save_path):
        print(f"Caminho de saída '{save_path}' não existe. Criando diretório...")
        os.makedirs(save_path)

    log_path = os.path.join(save_path, "log_erros.txt")
    with open(log_path, 'a', encoding='utf-8') as log_file:
        arquivos = [f for f in os.listdir(pasta) if f.endswith('.csv')]
        total_csvs = len(arquivos)

        print(f"\nPasta: {pasta}")
        print(f"Total de arquivos .csv encontrados: {total_csvs}")
        if total_csvs == 0:
            mensagem = f"Nenhum arquivo .csv encontrado na pasta: {pasta}\n"
            print(mensagem)
            log_file.write(mensagem)
            return

        cabecalho_referencia = None
        arquivos_validos_processados = 0
        nome_pasta = os.path.basename(os.path.abspath(pasta))
        nome_saida = os.path.join(save_path, f"{nome_pasta}_agrupado.csv")

        # Garante que um arquivo de uma execução anterior seja removido
        if os.path.exists(nome_saida):
            os.remove(nome_saida)
            print(f"Arquivo antigo '{nome_saida}' removido.")

        for arquivo in tqdm(arquivos, desc=f"Processando {nome_pasta}", leave=False):

            caminho_arquivo = os.path.join(pasta, arquivo)

            try:

                with open(caminho_arquivo, 'r', encoding='utf-8') as f:
                    linhas = f.readlines()

            except Exception as e:

                mensagem = f"Erro ao ler '{arquivo}': {e}\n"
                print(mensagem.strip())
                log_file.write(mensagem)
                continue

            if len(linhas) != 2:

                mensagem = f"Aviso: '{arquivo}' possui {len(linhas)} linhas. Ignorado.\n"
                print(mensagem.strip())
                log_file.write(mensagem)
                continue

            try:
                
                df = pd.read_csv(caminho_arquivo)

                # Se for o primeiro arquivo válido
                if cabecalho_referencia is None:
                    cabecalho_referencia = list(df.columns)
                    
                    df.to_csv(nome_saida, index=False, mode='w', header=True)
                    arquivos_validos_processados += 1
                
                # Para os arquivos seguintes
                else:

                    if list(df.columns) != cabecalho_referencia:

                        mensagem = f"Erro: Cabeçalho diferente em '{arquivo}'. Ignorado.\n"
                        print(mensagem.strip())
                        log_file.write(mensagem)

                        continue
                    
                    df.to_csv(nome_saida, index=False, mode='a', header=False)
                    arquivos_validos_processados += 1

            except Exception as e:

                mensagem = f"Erro ao processar ou salvar para '{arquivo}': {e}\n"
                print(mensagem.strip())
                log_file.write(mensagem)
        
        if arquivos_validos_processados > 0:

            print(f"Processo concluído. {arquivos_validos_processados} arquivos foram agrupados.")
            print(f"CSV final salvo como: {nome_saida}")

        else:

            mensagem = "Nenhum arquivo válido foi agrupado.\n"
            print(mensagem.strip())
            log_file.write(mensagem)

# --- Bloco de execução (sem alterações) ---
pastas_simples = [
    # "TESTE",
    # "fickett_score",
    "complex_networks"
    # "k-mer"
    # "orf"
]

for pasta in pastas_simples:
    agrupar_csvs(f'ZMays/{pasta}', 'ZMays/Agrup')