"""Executa ABCpred e BepiPred para uma proteína e compara as posições."""
import argparse
from datetime import datetime
import json
from pathlib import Path
import re
import sys
import shutil
import pandas as pd

from abcpred_scrapper import buscar_epitepos_abcpred, threshold, window, validar_parametros
from bepipred import buscar_epitopos_bepipred, preparar_sequencia
from relatorio import salvar_comparacao


def recomparar(origem, saida, limiar=None, score=None):
    origem = origem.resolve()
    fonte = json.loads((origem/'execucao.json').read_text(encoding='utf-8'))
    if fonte.get('estado', {}).get('abc') != 'concluido' or fonte.get('estado', {}).get('bepi') != 'concluido':
        raise ValueError('A pasta deve conter resultados concluídos dos dois preditores.')
    seq = preparar_sequencia(fonte['sequencia'])
    nome = fonte['proteina']
    if not re.fullmatch('[A-Za-z0-9_.-]+', nome):
        raise ValueError('Nome inválido no registro da execução.')
    limiar = fonte['bepi'].get('limiar_comparacao', .1512) if limiar is None else limiar
    score = fonte['bepi'].get('score_comparacao', 'bruto') if score is None else score
    abc = pd.read_csv(origem/'abc_resultado.csv')
    bepi = pd.read_csv(origem/'bepi_resultado.csv')
    from comparacao import comparar
    comparar(seq, abc, bepi, limiar, score)  # Validar antes de criar a pasta.
    pasta = saida / f'{nome}_recomparacao_{datetime.now().strftime("%Y%m%d_%H%M%S_%f")}'
    pasta.mkdir(parents=True)
    registro = {**fonte, 'bepi': {**fonte['bepi'], 'limiar_comparacao':limiar, 'score_comparacao':score},
                'modo':'recomparacao_local', 'origem':str(origem),
                'estado':{'abc':'concluido', 'bepi':'concluido', 'comparacao':'em_andamento'}, 'erros':{}}
    registro.pop('resumo', None)
    for arquivo in ['abc_resultado.csv', 'bepi_resultado.csv']:
        shutil.copy2(origem/arquivo, pasta/arquivo)
    (pasta/'entrada.fasta').write_text(f'>{nome}\n{seq}\n', encoding='utf-8')
    try:
        registro['resumo'] = salvar_comparacao(pasta, nome, seq, abc, bepi, limiar, score)
        registro['estado']['comparacao'] = 'concluida'
    except Exception as erro:
        registro['estado']['comparacao'] = 'falhou'
        registro['erros']['comparacao'] = str(erro)
        raise
    finally:
        (pasta/'execucao.json').write_text(json.dumps(registro, ensure_ascii=False, indent=2), encoding='utf-8')
    print(f'Recomparação local, sem envio aos serviços. Arquivos: {pasta.resolve()}')
    return pasta


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--nome', help='Nome da proteína (não busca sequências em bancos de dados)')
    e = p.add_mutually_exclusive_group()
    e.add_argument('--sequencia', help='Sequência de aminoácidos em uma linha')
    e.add_argument('--fasta', type=Path, help='Arquivo FASTA com uma única proteína')
    p.add_argument('--saida', type=Path, default=Path('resultados'), help='Pasta para as execuções')
    p.add_argument('--recomparar', type=Path, help='Pasta de execução existente; não consulta servidores')
    p.add_argument('--limiar-bepi', type=float, default=None)
    p.add_argument('--score-bepi', choices=['bruto', 'linear'], default=None)
    p.add_argument('--limiar-abc', type=float, default=threshold)
    p.add_argument('--janela-abc', type=int, choices=[10, 12, 14, 16, 18, 20], default=window)
    p.add_argument('--filtro-abc', choices=['off', 'on'], default='off')
    p.add_argument('--job-id-bepi', help='Retomar tarefa BepiPred já enviada, com a mesma sequência')
    p.add_argument('--tempo-limite', type=float, default=900, help='Espera máxima do BepiPred em segundos')
    a = p.parse_args()
    try:
        if a.recomparar:
            if a.nome or a.sequencia or a.fasta or a.job_id_bepi or any(
                t.split('=')[0] in ('--limiar-abc','--janela-abc','--filtro-abc','--tempo-limite') for t in sys.argv[1:]):
                raise ValueError('--recomparar aceita apenas --saida, --limiar-bepi e --score-bepi. '
                                 'Os peptídeos e parâmetros ABCpred originais são preservados.')
            recomparar(a.recomparar, a.saida, a.limiar_bepi, a.score_bepi)
            return 0
        a.limiar_bepi = .1512 if a.limiar_bepi is None else a.limiar_bepi
        a.score_bepi = 'bruto' if a.score_bepi is None else a.score_bepi
        nome = a.nome or input('Nome da proteína (ex.: Tp15): ').strip()
        if not re.fullmatch('[A-Za-z0-9_.-]+', nome):
            raise ValueError('Nome inválido. Use letras, números, ponto, hífen ou sublinhado.')
        texto = a.fasta.read_text(encoding='utf-8-sig') if a.fasta else a.sequencia or input('Sequência de aminoácidos: ')
        seq = preparar_sequencia(texto)
        validar_parametros(a.limiar_abc, a.janela_abc, a.filtro_abc)
        if len(seq) < a.janela_abc:
            raise ValueError('A sequência é menor que a janela ABCpred.')
        if not 0 <= a.limiar_bepi <= 1 or a.tempo_limite <= 0:
            raise ValueError('Limiar deve estar entre 0 e 1 e tempo limite deve ser positivo.')
        pasta = a.saida / f'{nome}_{datetime.now().strftime("%Y%m%d_%H%M%S_%f")}'
        pasta.mkdir(parents=True)
        (pasta / 'entrada.fasta').write_text(f'>{nome}\n{seq}\n', encoding='utf-8')
        registro = {'proteina': nome, 'sequencia': seq, 'abc': {'limiar': a.limiar_abc, 'janela': a.janela_abc, 'filtro': a.filtro_abc},
                    'bepi': {'versao': '3.0', 'limiar_comparacao': a.limiar_bepi, 'score_comparacao': a.score_bepi,
                             'limiar_servidor': 0.1512, 'top_epi': 0.2, 'roll_mean': 'no'}, 'estado': {}, 'erros': {}}
        dados = {}
        for sistema, executar in [('abc', lambda: buscar_epitepos_abcpred(nome, seq, limiar=a.limiar_abc,
                                                                       janela=a.janela_abc, filtro=a.filtro_abc)),
                                  ('bepi', lambda: buscar_epitopos_bepipred(nome, seq, job_id=a.job_id_bepi,
                                                                         limite_segundos=a.tempo_limite))]:
            try:
                df = executar()
                if df is None:
                    raise RuntimeError('Não foi possível interpretar o resultado recebido.')
                df.to_csv(pasta / f'{sistema}_resultado.csv', index=False, encoding='utf-8-sig')
                dados[sistema] = df
                registro['estado'][sistema] = 'concluido'
                if sistema == 'bepi':
                    registro['bepi'].update(df.attrs)
            except Exception as erro:
                registro['estado'][sistema] = 'falhou'
                registro['erros'][sistema] = str(erro)
                print(f'Falha em {sistema}: {erro}', file=sys.stderr)
            (pasta / 'execucao.json').write_text(json.dumps(registro, ensure_ascii=False, indent=2), encoding='utf-8')
        if len(dados) == 2:
            try:
                resumo = salvar_comparacao(pasta, nome, seq, dados['abc'], dados['bepi'], a.limiar_bepi, a.score_bepi)
                registro['resumo'] = resumo
                registro['estado']['comparacao'] = 'concluida'
            except Exception as erro:
                registro['estado']['comparacao'] = 'falhou'
                registro['erros']['comparacao'] = str(erro)
                print(f'Falha na comparação: {erro}', file=sys.stderr)
        else:
            registro['estado']['comparacao'] = 'nao_realizada'
        (pasta / 'execucao.json').write_text(json.dumps(registro, ensure_ascii=False, indent=2), encoding='utf-8')
        print(f'Arquivos salvos em: {pasta.resolve()}')
        return 1 if registro['erros'] else 0
    except (ValueError, OSError, KeyError, RuntimeError) as erro:
        print(f'Erro: {erro}', file=sys.stderr)
        return 1
    except KeyboardInterrupt:
        print('\nExecução interrompida. Resultados já salvos permanecem na pasta da execução.', file=sys.stderr)
        return 130


if __name__ == '__main__':
    raise SystemExit(main())
