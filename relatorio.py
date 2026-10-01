"""Salvamento compartilhado por predições novas e comparações locais."""
from pathlib import Path
from comparacao import comparar
from graficos import salvar_grafico


def salvar_comparacao(pasta, nome, seq, abc, bepi, limiar, score):
    pasta = Path(pasta)
    residuos, peptideos, regioes, resumo = comparar(seq, abc, bepi, limiar, score)
    for arquivo, df in [('comparacao_residuos', residuos), ('comparacao_peptideos', peptideos),
                         ('regioes_consenso', regioes)]:
        df.to_csv(pasta / f'{arquivo}.csv', index=False, encoding='utf-8-sig')
    salvar_grafico(nome, residuos, resumo, pasta/'grafico.html')
    print(f"Consenso: {resumo['residuos_consenso']}/{len(seq)} resíduos "
          f"({resumo['fracao_consenso']:.1%}), {resumo['regioes_consenso']} regiões.")
    return resumo
