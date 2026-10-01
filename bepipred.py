"""Envia uma proteína ao BepiPred 3.0 e salva scores por aminoácido."""
import argparse
import html
import re
import time
from io import BytesIO
from pathlib import Path
from urllib.parse import urljoin
from zipfile import ZipFile, is_zipfile

import pandas as pd
import requests

BASE = 'https://services.healthtech.dtu.dk'
CGI = BASE + '/cgi-bin/webface2.cgi'


def preparar_sequencia(texto):
    linhas = texto.strip().splitlines()
    if sum(linha.startswith('>') for linha in linhas) > 1:
        raise ValueError('Informe apenas uma proteína por execução.')
    if linhas and linhas[0].startswith('>'):
        linhas = linhas[1:]
    sequencia = ''.join(''.join(linhas).split()).upper()
    if not re.fullmatch('[ACDEFGHIKLMNPQRSTVWY]+', sequencia):
        raise ValueError('A sequência deve conter apenas aminoácidos canônicos, sem números.')
    if not 10 <= len(sequencia) <= 1023:
        raise ValueError('Use 10 a 1023 aminoácidos; o servidor trunca sequências maiores.')
    return sequencia


def ler_resultado(conteudo, sequencia):
    if not is_zipfile(BytesIO(conteudo)):
        raise ValueError('O conteúdo recebido não é um ZIP válido.')
    with ZipFile(BytesIO(conteudo)) as arquivo:
        nomes = [n for n in arquivo.namelist() if n.split('/')[-1] == 'raw_output.csv']
        if len(nomes) != 1:
            raise ValueError('raw_output.csv ausente ou ambíguo no resultado.')
        with arquivo.open(nomes[0]) as csv:
            df = pd.read_csv(csv, skipinitialspace=True)
    colunas = ['Residue', 'BepiPred-3.0 score', 'BepiPred-3.0 linear epitope score']
    if not set(colunas).issubset(df.columns):
        raise ValueError(f'Colunas inesperadas no resultado: {list(df.columns)}')
    if len(df) != len(sequencia) or ''.join(df['Residue'].astype(str)).upper() != sequencia:
        raise ValueError('Resultado incompleto ou sequência diferente da enviada.')
    df = df.rename(columns=dict(zip(colunas, ['Aminoácido', 'Score_Bruto_BP3', 'Score_Linear_BP3'])))
    for coluna in ['Score_Bruto_BP3', 'Score_Linear_BP3']:
        df[coluna] = pd.to_numeric(df[coluna], errors='raise')
        if not df[coluna].between(0, 1).all():
            raise ValueError(f'Scores inválidos em {coluna}.')
    df['Posição de Início'] = range(1, len(df) + 1)
    return df[['Posição de Início', 'Aminoácido', 'Score_Bruto_BP3', 'Score_Linear_BP3']]


def buscar_epitopos_bepipred(nome_proteina, sequencia_fasta, *, job_id=None,
                           limite_segundos=900, intervalo=10):
    """Retorna um DataFrame. job_id permite retomar sem reenviar a proteína.

    Interrupções e erros geram exceções; a submissão POST nunca é repetida
    automaticamente, pois pode ter sido aceita antes da queda de conexão.
    """
    sequencia = preparar_sequencia(sequencia_fasta)
    if not re.fullmatch(r'[A-Za-z0-9_.-]+', nome_proteina):
        raise ValueError('Use um nome sem espaços: letras, números, ponto, hífen ou sublinhado.')
    if job_id and not re.fullmatch(r'[A-Za-z0-9_-]+', job_id):
        raise ValueError('Identificador de tarefa inválido.')
    if limite_segundos <= 0 or intervalo <= 0:
        raise ValueError('Tempo limite e intervalo devem ser positivos.')
    with requests.Session() as sessao:
        sessao.headers.update({'User-Agent': 'Mozilla/5.0 (BepiPred Python client)'})
        if not job_id:
            payload = {'configfile': '/var/www/services/services/BepiPred-3.0/webface.cf',
                       'fasta': f'>{nome_proteina}\r\n{sequencia}\r\n', 'top_epi': '0.2',
                       'thr_epi': '0.1512', 'roll_mean': 'no'}
            print(f'Enviando {nome_proteina} ({len(sequencia)} aminoácidos) ao BepiPred...')
            resposta = sessao.post(CGI, data=payload,
                                  files={'uploadfile': ('', b'', 'application/octet-stream')}, timeout=60)
            resposta.raise_for_status()
            contexto = html.unescape(resposta.url + '\n' + resposta.text)
            match = re.search(r'(?:jobid\s*=\s*[\"\']?|/tmp/)([A-Za-z0-9_-]+)', contexto, re.I)
            if not match:
                raise RuntimeError('Servidor não informou uma tarefa. Resposta: ' + resposta.text[:300])
            job_id = match.group(1)
        pagina_url = CGI + f'?jobid={job_id}'
        print(f'Tarefa: {job_id}\nAcompanhamento: {pagina_url}', flush=True)
        prazo = time.monotonic() + limite_segundos
        ultimo_estado = None
        while time.monotonic() < prazo:
            try:
                restante = max(1, min(60, prazo - time.monotonic()))
                pagina = sessao.get(pagina_url, timeout=restante)
                pagina.raise_for_status()
                texto = html.unescape(pagina.text)
                link = re.search(r'href\s*=\s*[\"\']([^\"\']*bepipred3_results\.zip(?:\?[^\"\']*)?)[\"\']', texto, re.I)
                if link:
                    url_zip = urljoin(pagina.url, link.group(1))
                    resultado = sessao.get(url_zip, timeout=max(1, min(60, prazo-time.monotonic())))
                    if resultado.status_code in (202, 404):
                        time.sleep(min(intervalo, max(0, prazo-time.monotonic())))
                        continue
                    resultado.raise_for_status()
                    df = ler_resultado(resultado.content, sequencia)
                    df.attrs.update({'job_id': job_id, 'pagina_url': pagina_url, 'nome_proteina': nome_proteina})
                    print(f'Resultado recebido: {len(df)} aminoácidos.')
                    return df
                if 'ERROR: Could not find output file' in texto:
                    raise RuntimeError(f'O servidor não disponibilizou o resultado da tarefa {job_id}. Consulte {pagina_url}')
                match_estado = re.search(r"launchcheck\(['\"]([^'\"]+)", texto)
                estado = match_estado.group(1) if match_estado else 'aguardando página de resultado'
                if estado != ultimo_estado:
                    print('Estado: ' + {'queued': 'na fila', 'active': 'processando',
                                       'finished': 'processamento encerrado'}.get(estado, str(estado)), flush=True)
                    ultimo_estado = estado
                if estado in ('failed', 'error', 'killed', 'rejected', 'expired', 'deleted'):
                    raise RuntimeError(f'Tarefa {job_id} encerrada com estado {estado}. Consulte {pagina_url}')
            except (requests.ConnectionError, requests.Timeout) as erro:
                print(f'Conexão interrompida ({type(erro).__name__}); tentando consultar a mesma tarefa novamente.', flush=True)
            except requests.HTTPError as erro:
                if erro.response is None or erro.response.status_code not in (429, 500, 502, 503, 504):
                    raise
                print('Servidor temporariamente indisponível; aguardando para consultar a mesma tarefa.', flush=True)
            time.sleep(min(intervalo, max(0, prazo - time.monotonic())))
        raise TimeoutError(f'Tempo limite atingido. Retome usando --job-id {job_id}. Página: {pagina_url}')


# Compatibilidade com o nome usado pelo código anterior.
buscar_epitepos_bepipred = buscar_epitopos_bepipred


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--nome', help='Nome da proteína, por exemplo Tp15')
    entrada = parser.add_mutually_exclusive_group()
    entrada.add_argument('--sequencia', help='Sequência de aminoácidos')
    entrada.add_argument('--fasta', type=Path, help='Arquivo FASTA com uma proteína')
    parser.add_argument('--saida', type=Path, help='Caminho do CSV de saída')
    parser.add_argument('--job-id', help='Retomar tarefa já enviada')
    parser.add_argument('--tempo-limite', type=float, default=900, help='Espera máxima em segundos (padrão 900)')
    args = parser.parse_args()
    try:
        nome = args.nome or input('Nome da proteína (ex.: Tp15): ').strip()
        if args.fasta:
            texto = args.fasta.read_text(encoding='utf-8-sig')
        else:
            texto = args.sequencia or input('Sequência de aminoácidos (em uma linha): ').strip()
        df = buscar_epitopos_bepipred(nome, texto, job_id=args.job_id, limite_segundos=args.tempo_limite)
        caminho = args.saida or Path('resultados') / f'{nome}_bepipred.csv'
        caminho.parent.mkdir(parents=True, exist_ok=True)
        df.to_csv(caminho, index=False, encoding='utf-8-sig')
        print(f'CSV salvo em: {caminho.resolve()}')
        print(df.head(10).to_string(index=False))
    except (ValueError, RuntimeError, TimeoutError, requests.RequestException, OSError) as erro:
        print(f'Erro: {erro}')
        return 1
    except KeyboardInterrupt:
        print('\nConsulta interrompida. A tarefa enviada pode continuar no servidor; guarde o identificador acima.')
        return 130
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
