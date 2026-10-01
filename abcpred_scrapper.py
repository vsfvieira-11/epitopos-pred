import requests
import pandas as pd
from io import StringIO

threshold = 0.51
window = 16
filter = 0 #0 = off, 1 = on

def buscar_epitepos_abcpred(nome_proteina, sequencia_fasta):
    # ATENÇÃO: Confirme se essa é a URL exata do 'action' do formulário inspecionando a página (F12)
    url_action = "https://webs.iiitd.edu.in/cgibin/abcpred/test1_main.pl" 
    
    payload = {
        "SEQNAME": nome_proteina,
        "SEQ": sequencia_fasta,
        "Threshold": str(threshold),
        "window": str(window),
        "filter": "off" if filter == 0 else "on",
        "submit": "Submit sequence"
    }

    # Simulando um navegador real (alguns servidores bloqueiam requests sem User-Agent)
    headers = {
        'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/120.0.0.0 Safari/537.36'
    }


    print("Enviando dados para o ABCpred...")
    resposta = requests.post(url_action, data=payload, headers=headers)
    
    try:
        # Pede pro Pandas extrair TODAS as tabelas do HTML retornado
        # Usamos StringIO para converter o texto da resposta em um formato que o Pandas prefere
        tabelas = pd.read_html(StringIO(resposta.text))
        
        # Como pode haver mais de uma tabela (a de Input Info e a de Resultados),
        # vamos varrer a lista e pegar a que tem a coluna "Rank" ou "Score"
        df_resultados = None
        for df in tabelas:
            # Pandas converte a primeira linha em cabeçalho, vamos checar se é a tabela certa
            if 'Rank' in df.columns or 'Score' in df.columns or (0 in df.columns and 'Rank' in str(df[0].values)):
                # Se não vier com cabeçalho limpo, forçamos a primeira linha a ser o cabeçalho
                if 'Rank' not in df.columns:
                    df.columns = df.iloc[0] # Define a 1ª linha como nome das colunas
                    df = df[1:] # Remove a 1ª linha dos dados
                    df = df.reset_index(drop=True)
                
                df_resultados = df
                break
        
        # Se mesmo assim não achou por nome, tenta pegar a última tabela da página (fallback)
        if df_resultados is None and len(tabelas) > 0:
            df_resultados = tabelas[-1]

        # Limpeza básica: garante que a pontuação é número (float) e a posição é inteiro (int)
        if df_resultados is not None:
            # Pega explicitamente APENAS as 4 primeiras colunas (ignora a 5ª fantasma)
            df_resultados = df_resultados.iloc[:, :4]
            
            # Agora sim renomeamos com segurança (pois garantimos que há 4 colunas)
            df_resultados.columns = ['Rank', 'Sequência', 'Posição de Início', 'Score']
            
            # Converte os tipos numéricos
            df_resultados['Score'] = pd.to_numeric(df_resultados['Score'], errors='coerce')
            df_resultados['Posição de Início'] = pd.to_numeric(df_resultados['Posição de Início'], errors='coerce')
            
            # Remove as linhas nulas E reseta o index na mesma linha!
            df_resultados = df_resultados.dropna(subset=['Posição de Início']).reset_index(drop=True)
            
            # Converte a posição inicial para inteiro
            df_resultados['Posição de Início'] = df_resultados['Posição de Início'].astype(int)
            
        return df_resultados

    except Exception as e:
        print(f"Erro ao tentar extrair a tabela HTML: {e}")
        return None

# --- ÁREA DE TESTE ---
if __name__ == "__main__":
    seq_teste = "MVKRGRFALCLAVLLGACSFSSIPNGTYRATYQDFDENGWKDFLEVTFDGGKMVQVVYDYQHKEGRFKSQDADYHRVMYASSGIGPEKAFRELADALLEKGNPEMVDVVTGATVSSQSFRRLGRALLQSARRGEKEAIISR"
    
    # Chama a função, que agora nos retorna um DataFrame prontinho!
    df_abcpred = buscar_epitepos_abcpred("Tp15_Teste", seq_teste)
    
    if df_abcpred is not None and not df_abcpred.empty:
        print("\nSUCESSO TOTAL! A Geovana vai pirar. Olha os dados aí:")
        print(df_abcpred) # Mostra as 10 primeiras linhas da tabela
        
        # Como é um DataFrame, você já pode salvar em CSV se quiser:
        # df_abcpred.to_csv("resultados_abcpred.csv", index=False)
    else:
        print("\nO Pandas não conseguiu encontrar nenhuma tabela no HTML retornado.")