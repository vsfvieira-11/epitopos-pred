import requests
import re
import time
import pandas as pd
from io import BytesIO
from zipfile import ZipFile

def buscar_epitepos_bepipred(nome_proteina, sequencia_fasta):
    url_cgi = "https://services.healthtech.dtu.dk/cgi-bin/webface2.cgi"
    
    payload = {
        "configfile": "/var/www/services/services/BepiPred-3.0/webface.cf",
        "fasta": f">{nome_proteina}\n{sequencia_fasta}",
        "top_epi": "0.2",
        "thr_epi": "0.1512",
        "roll_mean": "yes"
    }

    headers = {
        'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36'
    }

    print("Enviando a sequência para o servidor ultrarrápido da DTU...")
    sessao = requests.Session()
    resposta = sessao.post(url_cgi, data=payload, headers=headers)

    match = re.search(r'jobid=([A-Za-z0-9_]+)', resposta.text, re.IGNORECASE)
    if not match:
        print("Erro: Não foi possível capturar o Job ID no servidor.")
        return None
        
    job_id = match.group(1)
    print(f"Job ID capturado: {job_id} - Aguardando processamento...")
    
    # A URL direta do ZIP (A grande sacada que você pescou no HTML)
    url_zip = f"https://services.healthtech.dtu.dk/services/BepiPred-3.0/tmp/{job_id}/bepipred3_results.zip"
    
    tentativas = 0
    max_tentativas = 40 # Cerca de 3 minutos de limite
    
    # Em vez de olhar o HTML, ficamos consultando o link do ZIP diretamente.
    while tentativas < max_tentativas:
        print(f"A processar na GPU da Dinamarca... (Tentativa {tentativas+1}/{max_tentativas})", end="\r")
        
        # Tenta aceder ao ficheiro
        res_zip = sessao.get(url_zip, headers=headers)
        
        if res_zip.status_code == 200:
            print("\n\nProcesso concluído! O arquivo ZIP foi gerado. A extrair CSV...")
            
            # 1. Carrega o ZIP para a RAM e abre-o
            with ZipFile(BytesIO(res_zip.content)) as arquivo_zip:
                # 2. Abre só o CSV que nos interessa
                with arquivo_zip.open('raw_output.csv') as arquivo_csv:
                    df_bepipred = pd.read_csv(arquivo_csv)
            
            # Limpeza e formatação para ficar amigável
            df_bepipred = df_bepipred.rename(columns={
                'Residue': 'Aminoácido',
                'BepiPred-3.0 score': 'Score_Bruto_BP3',
                'BepiPred-3.0 linear epitope score': 'Score_Consenso_BP3'
            })
            
            # Recria a numeração (o index começa em 0, somamos 1)
            df_bepipred['Posição de Início'] = df_bepipred.index + 1
            
            # Devolve exatamente as colunas prontas para o cruzamento
            return df_bepipred[['Posição de Início', 'Aminoácido', 'Score_Bruto_BP3', 'Score_Consenso_BP3']]
            
        elif res_zip.status_code == 404:
            # O ZIP não existe ainda, o modelo está a correr.
            time.sleep(5)
            tentativas += 1
        else:
            print(f"\nErro inesperado no servidor. Código HTTP: {res_zip.status_code}")
            return None
            
    print("\nTempo limite excedido. O servidor demorou muito para gerar o arquivo.")
    return None

# --- ÁREA DE TESTE ---
if __name__ == "__main__":
    seq_teste = "MVKRGRFALCLAVLLGACSFSSIPNGTYRATYQDFDENGWKDFLEVTFDGGKMVQVVYDYQHKEGRFKSQDADYHRVMYASSGIGPEKAFRELADALLEKGNPEMVDVVTGATVSSQSFRRLGRALLQSARRGEKEAIISR"
    
    df_bp = buscar_epitepos_bepipred("Tp15_Teste", seq_teste)
    
    if df_bp is not None:
        print("\nSUCESSO TOTAL! A tabela CSV foi extraída e está pronta:")
        print(df_bp.head(10))