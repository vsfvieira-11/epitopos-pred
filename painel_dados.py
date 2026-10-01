"""Operações locais do painel, independentes da interface Streamlit."""
from datetime import datetime
from io import BytesIO
import json
import math
from pathlib import Path
from zipfile import ZipFile, ZIP_DEFLATED
import pandas as pd
from comparacao import comparar


def listar_execucoes(raiz):
    validas, problemas = [], []
    raiz = Path(raiz)
    if not raiz.is_dir():
        return [], [f'Pasta não encontrada: {raiz}']
    for pasta in sorted(raiz.iterdir(), reverse=True):
        if not pasta.is_dir() or not (pasta/'execucao.json').exists():
            continue
        try:
            registro = json.loads((pasta/'execucao.json').read_text(encoding='utf-8'))
            if not isinstance(registro, dict) or not registro.get('proteina') or not registro.get('sequencia'):
                raise ValueError('Registro sem proteína ou sequência')
            if not isinstance(registro.get('abc'), dict) or not isinstance(registro.get('bepi'), dict):
                raise ValueError('Parâmetros dos preditores ausentes')
            limiar = float(registro['abc']['limiar'])
            limiar_bp = float(registro['bepi'].get('limiar_comparacao',.1512))
            if not math.isfinite(limiar) or not .1 <= limiar <= 1 or not math.isfinite(limiar_bp) or not 0 <= limiar_bp <= 1:
                raise ValueError('Limiar inválido no registro')
            if 'janela' not in registro['abc'] or 'filtro' not in registro['abc']:
                raise ValueError('Janela ou filtro ABCpred ausente')
            if any(registro.get('estado', {}).get(p) != 'concluido' for p in ['abc', 'bepi']):
                raise ValueError('Um dos preditores não foi concluído')
            if any(not (pasta/f'{p}_resultado.csv').is_file() for p in ['abc','bepi']):
                raise ValueError('CSV de predição ausente')
            validas.append({'pasta':pasta, 'registro':registro})
        except (ValueError, OSError, TypeError, KeyError, AttributeError) as erro:
            problemas.append(f'{pasta.name}: {erro}')
    return validas, problemas


def carregar_execucao(pasta):
    pasta = Path(pasta)
    registro = json.loads((pasta/'execucao.json').read_text(encoding='utf-8'))
    abc = pd.read_csv(pasta/'abc_resultado.csv')
    bepi = pd.read_csv(pasta/'bepi_resultado.csv')
    # Validar todos os peptídeos antes de permitir filtros que poderiam ocultar erros.
    comparar(registro['sequencia'], abc, bepi)
    return registro, abc, bepi


def analisar_local(registro, abc, bepi, score, limiar_bp, limiar_abc):
    limiar_origem = float(registro['abc']['limiar'])
    if not limiar_origem <= limiar_abc <= 1:
        raise ValueError('O corte ABCpred local não pode ser menor que o limiar enviado ao servidor.')
    selecionados = abc.loc[pd.to_numeric(abc['Score'], errors='raise') >= limiar_abc].copy()
    residuos, peptideos, regioes, resumo = comparar(registro['sequencia'], selecionados, bepi, limiar_bp, score)
    peptideos['Tamanho'] = peptideos['Fim'] - peptideos['Inicio'] + 1
    resumo['limiar_abc_local'] = limiar_abc
    resumo['peptideos_abc_originais'] = len(abc)
    return residuos, peptideos, regioes, resumo


def filtrar_candidatos(peptideos, inicio=1, fim=None, fracao_min=0, tamanho_min=1, score_min=0, busca=''):
    if fim is None:
        fim = int(peptideos['Fim'].max()) if not peptideos.empty else inicio
    mask = ((peptideos['Fim'] >= inicio) & (peptideos['Inicio'] <= fim) &
            (peptideos['Fracao_BP_positiva'] >= fracao_min) &
            (peptideos['Tamanho'] >= tamanho_min) & (peptideos['Score_ABC'] >= score_min))
    if busca.strip():
        mask &= peptideos['Sequencia'].str.contains(busca.strip().upper(), regex=False)
    return peptideos.loc[mask].copy().reset_index(drop=True)


def csv_bytes(df):
    return df.to_csv(index=False).encode('utf-8-sig')


def pacote_analise(registro, origem, residuos, peptideos, regioes, resumo, filtros, candidatos):
    """Exportação analítica; parâmetros de análise e de visualização ficam separados."""
    buffer = BytesIO()
    with ZipFile(buffer, 'w', ZIP_DEFLATED) as arquivo:
        for nome, df in [('comparacao_residuos',residuos), ('comparacao_peptideos',peptideos),
                         ('regioes_consenso',regioes), ('candidatos_filtrados',candidatos)]:
            arquivo.writestr(nome+'.csv', csv_bytes(df))
        arquivo.writestr('entrada.fasta', f">{registro['proteina']}\n{registro['sequencia']}\n")
        arquivo.writestr('parametros.json', json.dumps({'proteina':registro['proteina'],
            'gerado_em':datetime.now().isoformat(), 'origem':str(origem),
            'predicao_original':{'abc':registro['abc'],'bepi':registro['bepi']},
            'analise_atual':resumo, 'filtros_de_visualizacao':filtros,
            'nota':'Filtros de candidatos não recalculam o consenso. Análise exploratória de epítopos preditos.'},
            ensure_ascii=False, indent=2))
    return buffer.getvalue()
