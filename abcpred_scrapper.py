import requests
import pandas as pd
from io import StringIO
import argparse
import math
from pathlib import Path
import re

threshold = 0.51
window = 16
filter = 0 #0 = off, 1 = on

def validar_parametros(limiar=threshold, janela=window, filtro='off'):
    if not math.isfinite(limiar) or not 0.1 <= limiar <= 1:
        raise ValueError('Limiar ABCpred deve estar entre 0.1 e 1.0.')
    if janela not in (10, 12, 14, 16, 18, 20):
        raise ValueError('Janela ABCpred deve ser 10, 12, 14, 16, 18 ou 20.')
    if filtro not in ('off', 'on'):
        raise ValueError('Filtro ABCpred deve ser off ou on.')


def buscar_epitepos_abcpred(nome_proteina, sequencia_fasta, *, limiar=threshold,
                           janela=window, filtro='off', timeout=120):
    validar_parametros(limiar, janela, filtro)
    if not math.isfinite(timeout) or timeout <= 0:
        raise ValueError('Timeout deve ser positivo e finito.')
    sequencia_fasta = ''.join(sequencia_fasta.split()).upper()
    if not re.fullmatch('[ACDEFGHIKLMNPQRSTVWY]+', sequencia_fasta) or len(sequencia_fasta) < janela:
        raise ValueError('Informe aminoácidos canônicos e uma sequência com pelo menos o tamanho da janela.')
    # ATENÇÃO: Confirme se essa é a URL exata do 'action' do formulário inspecionando a página (F12)
    url_action = "https://webs.iiitd.edu.in/cgibin/abcpred/test1_main.pl" 
    
    payload = {
        "SEQNAME": nome_proteina,
        "SEQ": sequencia_fasta,
        "Threshold": str(limiar),
        "window": str(janela),
        "filter": filtro,
        "submit": "Submit sequence"
    }

    # Simulando um navegador real (alguns servidores bloqueiam requests sem User-Agent)
    headers = {
        'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/120.0.0.0 Safari/537.36'
    }


    print("Enviando dados para o ABCpred...")
    resposta = requests.post(url_action, data=payload, headers=headers, timeout=timeout)
    resposta.raise_for_status()
    
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

buscar_epitopos_abcpred = buscar_epitepos_abcpred


def main():
    p = argparse.ArgumentParser(description='Consulta ABCpred e salva seus peptídeos em CSV.')
    p.add_argument('--nome')
    e = p.add_mutually_exclusive_group()
    e.add_argument('--sequencia')
    e.add_argument('--fasta', type=Path)
    p.add_argument('--limiar', type=float, default=threshold)
    p.add_argument('--janela', type=int, choices=[10, 12, 14, 16, 18, 20], default=window)
    p.add_argument('--filtro', choices=['off', 'on'], default='off')
    p.add_argument('--timeout', type=float, default=120)
    p.add_argument('--saida', type=Path)
    a = p.parse_args()
    try:
        nome = a.nome or input('Nome da proteína: ').strip()
        if not re.fullmatch('[A-Za-z0-9_.-]+', nome):
            raise ValueError('Nome inválido; use letras, números, ponto, hífen ou sublinhado.')
        texto = a.fasta.read_text(encoding='utf-8-sig') if a.fasta else a.sequencia or input('Sequência de aminoácidos: ')
        linhas = texto.strip().splitlines()
        if sum(l.startswith('>') for l in linhas) > 1:
            raise ValueError('Informe apenas uma proteína.')
        if linhas and linhas[0].startswith('>'):
            linhas = linhas[1:]
        seq = ''.join(''.join(linhas).split()).upper()
        df = buscar_epitepos_abcpred(nome, seq, limiar=a.limiar, janela=a.janela,
                                   filtro=a.filtro, timeout=a.timeout)
        if df is None:
            raise RuntimeError('Não foi possível interpretar a resposta ABCpred.')
        destino = a.saida or Path('resultados') / f'{nome}_abcpred.csv'
        destino.parent.mkdir(parents=True, exist_ok=True)
        df.to_csv(destino, index=False, encoding='utf-8-sig')
        print(f'Parâmetros: limiar={a.limiar}, janela={a.janela}, filtro={a.filtro}')
        print(f'{len(df)} peptídeos. CSV salvo em: {destino.resolve()}')
        print(df.head(10).to_string(index=False))
    except (ValueError, RuntimeError, OSError, requests.RequestException) as erro:
        print(f'Erro: {erro}')
        return 1
    except KeyboardInterrupt:
        print('\nConsulta interrompida.')
        return 130
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
