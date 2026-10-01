"""Compara peptídeos ABCpred com scores BepiPred da mesma sequência."""
import math
import pandas as pd


def comparar(sequencia, abc, bepi, limiar=0.1512, score='bruto'):
    if not math.isfinite(limiar) or not 0 <= limiar <= 1:
        raise ValueError('Limiar BepiPred deve estar entre 0 e 1.')
    if score not in ('bruto', 'linear'):
        raise ValueError('Score BepiPred deve ser bruto ou linear.')
    if abc is None or bepi is None:
        raise ValueError('Resultados dos dois preditores são necessários.')
    bp = bepi.copy()
    col = 'Score_Bruto_BP3' if score == 'bruto' else 'Score_Linear_BP3'
    pos = pd.to_numeric(bp['Posição de Início'], errors='raise')
    if len(bp) != len(sequencia) or list(pos) != list(range(1, len(sequencia)+1)):
        raise ValueError('BepiPred não cobre as posições da sequência na ordem esperada.')
    if ''.join(bp['Aminoácido'].astype(str)).upper() != sequencia:
        raise ValueError('A sequência do BepiPred difere da sequência de entrada.')
    for c in ('Score_Bruto_BP3', 'Score_Linear_BP3'):
        bp[c] = pd.to_numeric(bp[c], errors='raise')
        if not bp[c].between(0, 1).all():
            raise ValueError('Scores BepiPred inválidos.')
    positivos = list(bp[col] >= limiar)
    cobertura = [0] * len(sequencia)
    max_scores = [None] * len(sequencia)
    peptideos = []
    for _, r in abc.iterrows():
        inicio_n = float(r['Posição de Início'])
        if not math.isfinite(inicio_n) or inicio_n != int(inicio_n):
            raise ValueError('Posição ABCpred não inteira.')
        inicio = int(inicio_n)
        peptideo = str(r['Sequência']).strip().upper()
        fim = inicio + len(peptideo) - 1
        if not peptideo or not 1 <= inicio <= fim <= len(sequencia):
            raise ValueError('Intervalo ABCpred fora da proteína.')
        if sequencia[inicio-1:fim] != peptideo:
            raise ValueError(f'Peptídeo ABCpred na posição {inicio} não corresponde à proteína.')
        abc_score = float(r['Score'])
        if not math.isfinite(abc_score) or not 0 <= abc_score <= 1:
            raise ValueError('Score ABCpred inválido.')
        for i in range(inicio-1, fim):
            cobertura[i] += 1
            max_scores[i] = abc_score if max_scores[i] is None else max(max_scores[i], abc_score)
        n = sum(positivos[inicio-1:fim])
        peptideos.append({'Rank_ABC': r['Rank'], 'Inicio': inicio, 'Fim': fim, 'Sequencia': peptideo,
                         'Score_ABC': abc_score, 'Score_BP_medio': float(bp[col].iloc[inicio-1:fim].mean()),
                         'Residuos_BP_positivos': n, 'Fracao_BP_positiva': n/len(peptideo)})
    residuos = pd.DataFrame({'Posicao': range(1, len(sequencia)+1), 'Aminoacido': list(sequencia),
                             'Score_BP_bruto': list(bp['Score_Bruto_BP3']),
                             'Score_BP_linear': list(bp['Score_Linear_BP3']),
                             'BP_positivo': positivos, 'ABC_coberto': [n > 0 for n in cobertura],
                             'ABC_n_peptideos': cobertura, 'ABC_score_max': max_scores,
                             'Consenso': [b and n > 0 for b, n in zip(positivos, cobertura)]})
    regioes = []
    inicio = None
    for i, positivo in enumerate(list(residuos['Consenso']) + [False], 1):
        if positivo and inicio is None:
            inicio = i
        elif not positivo and inicio is not None:
            fim = i-1
            regioes.append({'Inicio': inicio, 'Fim': fim, 'Tamanho': fim-inicio+1,
                            'Sequencia': sequencia[inicio-1:fim]})
            inicio = None
    ambos = int(residuos['Consenso'].sum())
    uniao = sum(b or n > 0 for b, n in zip(positivos, cobertura))
    resumo = {'tamanho_proteina': len(sequencia), 'peptideos_abc': len(peptideos),
              'score_bp_usado': score, 'limiar_bp': limiar, 'regra_bp': 'score >= limiar',
              'residuos_abc': sum(n > 0 for n in cobertura), 'residuos_bp': sum(positivos),
              'residuos_consenso': ambos, 'fracao_consenso': ambos/len(sequencia),
              'regioes_consenso': len(regioes), 'jaccard_cobertura': ambos/uniao if uniao else None}
    return (residuos, pd.DataFrame(peptideos, columns=['Rank_ABC', 'Inicio', 'Fim', 'Sequencia', 'Score_ABC',
                                                    'Score_BP_medio', 'Residuos_BP_positivos', 'Fracao_BP_positiva']),
            pd.DataFrame(regioes, columns=['Inicio', 'Fim', 'Tamanho', 'Sequencia']), resumo)
