"""Interface local: streamlit run app.py"""
import os
from pathlib import Path
import streamlit as st
from graficos import criar_figura
from painel_dados import (listar_execucoes, carregar_execucao, analisar_local,
                          filtrar_candidatos, csv_bytes, pacote_analise)

ROOT = Path(__file__).resolve().parent


def tabela(df, key=None, selecionar=False):
    config = {
        'Inicio':'Início', 'Fim':'Fim', 'Tamanho':'Comprimento (aa)', 'Sequencia':'Sequência',
        'Posicao':'Posição', 'Aminoacido':'Aminoácido', 'Rank_ABC':'Rank ABCpred',
        'Score_ABC':st.column_config.NumberColumn('Score ABCpred', format='%.3f'),
        'Score_BP_medio':st.column_config.NumberColumn('Média BepiPred', format='%.4f'),
        'Fracao_BP_positiva':st.column_config.ProgressColumn('Fração BepiPred positiva', min_value=0, max_value=1, format='%.2f'),
        'Residuos_BP_positivos':'Resíduos BepiPred positivos', 'Consenso':'Consenso',
        'ABC_coberto':'Coberto pelo ABCpred', 'BP_positivo':'BepiPred positivo',
        'ABC_n_peptideos':'Peptídeos ABCpred sobrepostos',
        'Score_BP_bruto':st.column_config.NumberColumn('BepiPred bruto',format='%.4f'),
        'Score_BP_linear':st.column_config.NumberColumn('BepiPred linear',format='%.4f'),
        'ABC_score_max':st.column_config.NumberColumn('Maior score ABCpred sobreposto',format='%.3f')}
    return st.dataframe(df, width='stretch', hide_index=True, column_config=config, key=key,
                        on_select='rerun' if selecionar else 'ignore', selection_mode='single-row')


def main():
    st.set_page_config(page_title='Epítopos · Painel de análise', page_icon='🧬', layout='wide')
    st.title('Explorador de epítopos')
    st.caption('ABCpred + BepiPred 3.0 · Comparação de candidatos a epítopos de células B')
    with st.sidebar:
        st.header('Dados da análise')
        raiz = st.text_input('Pasta de resultados', value=os.environ.get('EPITOPOS_RESULTADOS',str(ROOT/'resultados')))
        st.button('Atualizar execuções')
    execucoes, problemas = listar_execucoes(raiz)
    if not execucoes:
        st.info('Nenhuma execução completa encontrada. Gere resultados com o main.py ou selecione outra pasta.')
        st.code('python main.py --nome Tp15 --fasta exemplos/Tp15.fasta --score-bepi linear')
        if problemas:
            st.write(problemas)
        return
    with st.sidebar:
        nomes = sorted({e['registro']['proteina'] for e in execucoes})
        proteina = st.selectbox('Proteína', nomes)
        disponiveis = [e for e in execucoes if e['registro']['proteina']==proteina]
        pasta = st.selectbox('Execução', [e['pasta'] for e in disponiveis], format_func=lambda p:p.name)
    try:
        registro, abc, bepi = carregar_execucao(pasta)
    except (ValueError, KeyError, OSError, TypeError) as erro:
        st.error(f'Não foi possível carregar esta execução: {erro}')
        return
    chave = pasta.name
    with st.sidebar:
        st.subheader('Recalcular comparação')
        modo_original = registro['bepi'].get('score_comparacao','bruto')
        modo = st.radio('Score BepiPred', ['linear','bruto'],
                        index=0 if modo_original=='linear' else 1, key='modo_'+chave,
                        help='Linear suaviza scores na sequência. Bruto usa o score de cada resíduo.')
        limiar = st.number_input('Limiar BepiPred', min_value=0.0, max_value=1.0,
                                value=float(registro['bepi'].get('limiar_comparacao',.1512)),
                                step=.001, format='%.4f', key='bp_'+chave)
        limiar_origem = float(registro['abc']['limiar'])
        corte_abc = st.number_input('Corte ABCpred local', min_value=limiar_origem, max_value=1.0,
                                   value=limiar_origem, step=.01, format='%.2f', key='abc_'+chave,
                                   help='Remove peptídeos salvos com score abaixo do corte e recalcula cobertura e consenso. Não permite recuperar peptídeos não retornados pelo servidor.')
        st.caption('Estes controles recalculam os dados localmente. Janela e filtro ABCpred permanecem os da execução original.')
    try:
        residuos, peptideos, regioes, resumo = analisar_local(registro,abc,bepi,modo,limiar,corte_abc)
    except (ValueError, KeyError, TypeError) as erro:
        st.error(f'Não foi possível comparar os dados: {erro}')
        return
    st.subheader(proteina)
    st.caption(f"{pasta.name} · ABCpred: janela {registro['abc']['janela']}, filtro {registro['abc']['filtro']} · "
               f"Regra atual: BepiPred {modo} ≥ {limiar:g} e ABCpred ≥ {corte_abc:g}")
    a,b,c,d = st.columns(4)
    a.metric('Comprimento',f"{resumo['tamanho_proteina']} aa")
    b.metric('Peptídeos ABCpred',resumo['peptideos_abc'])
    b.caption(f'{len(abc)} peptídeos na execução original')
    c.metric('Resíduos em consenso',str(resumo['residuos_consenso']))
    c.caption(f"{resumo['fracao_consenso']:.1%} da proteína")
    d.metric('Regiões',resumo['regioes_consenso'])
    st.caption('Consenso = resíduo coberto pelo ABCpred e score BepiPred acima ou igual ao limiar. Não representa confirmação experimental.')
    visao, candidatos, aminoacidos, trechos, origem = st.tabs(['Visão geral','Candidatos','Aminoácidos','Regiões de consenso','Origem e exportação'])
    figura = criar_figura(proteina,residuos,resumo)
    with visao:
        st.plotly_chart(figura, width='stretch', key='grafico_'+chave,
                        config={'displaylogo':False,'toImageButtonOptions':{'format':'svg','filename':proteina+'_comparacao'}})
        st.caption('Arraste para ampliar, use a barra inferior para navegar e a câmera para exportar SVG. Ocultar uma curva não altera o consenso.')
    with candidatos:
        st.subheader('Peptídeos candidatos')
        st.write('Filtre os peptídeos e selecione uma linha para inspecionar seus aminoácidos e sua posição no gráfico.')
        x,y,z = st.columns(3)
        fracao = x.slider('Fração mínima de resíduos BepiPred positivos (%)',0,100,0,key='fracao_'+chave)
        tamanho = y.number_input('Comprimento mínimo do peptídeo (aa)',1,len(registro['sequencia']),1,key='tam_'+chave)
        busca = z.text_input('Buscar sequência no peptídeo',key='busca_'+chave)
        intervalo = st.slider('Intervalo de posições para procurar candidatos',1,len(registro['sequencia']),
                              (1,len(registro['sequencia'])),key='intervalo_'+chave)
        ordenacao = st.selectbox('Ordenar por',['Posição','Fração BepiPred positiva','Score ABCpred','Média BepiPred'],key='ord_'+chave)
        filtrados = filtrar_candidatos(peptideos,intervalo[0],intervalo[1],fracao/100,int(tamanho),busca=busca)
        coluna = {'Posição':'Inicio','Fração BepiPred positiva':'Fracao_BP_positiva','Score ABCpred':'Score_ABC','Média BepiPred':'Score_BP_medio'}[ordenacao]
        filtrados = filtrados.sort_values(coluna,ascending=ordenacao=='Posição',kind='stable').reset_index(drop=True)
        filtros = {'fracao_bp_minima':fracao/100,'comprimento_minimo':int(tamanho),'busca_sequencia':busca,
                   'intervalo':list(intervalo),'ordenacao':ordenacao}
        st.caption(f'{len(filtrados)} de {len(peptideos)} peptídeos exibidos. Intervalo seleciona peptídeos que o intersectam. '
                   'Os filtros desta tabela não alteram os totais nem as regiões de consenso.')
        evento = tabela(filtrados,key=f'candidatos_{chave}_{modo}_{limiar}_{corte_abc}_{filtros}',selecionar=True)
        st.download_button('Baixar candidatos filtrados (CSV)',csv_bytes(filtrados),f'{proteina}_candidatos.csv','text/csv')
        if filtrados.empty:
            st.info('Nenhum candidato atende aos filtros atuais. Experimente reduzir a cobertura mínima ou o comprimento.')
        elif evento.selection.rows:
            i = evento.selection.rows[0]
            if i < len(filtrados):
                r = filtrados.iloc[i]
                st.subheader(f"Peptídeo {int(r['Inicio'])}–{int(r['Fim'])}")
                st.code(r['Sequencia'],language=None)
                detalhe = residuos.loc[residuos['Posicao'].between(int(r['Inicio']),int(r['Fim']))]
                figura_detalhe = criar_figura(proteina,residuos,resumo)
                figura_detalhe.add_vrect(x0=int(r['Inicio'])-.5,x1=int(r['Fim'])+.5,
                                        fillcolor='#f59e0b',opacity=.15,line_width=0,row=1,col=1)
                figura_detalhe.update_xaxes(range=[int(r['Inicio'])-.5,int(r['Fim'])+.5])
                st.plotly_chart(figura_detalhe,width='stretch',key='detalhe_'+chave)
                tabela(detalhe)
    with aminoacidos:
        st.subheader('Tabela por aminoácido')
        inicio,fim = st.slider('Posições exibidas',1,len(registro['sequencia']),
                               (1,len(registro['sequencia'])),key='res_intervalo_'+chave)
        apenas = st.checkbox('Mostrar somente resíduos em consenso',key='consenso_'+chave)
        aa = st.multiselect('Aminoácidos', sorted(set(registro['sequencia'])),key='aa_'+chave)
        linhas = residuos.loc[residuos['Posicao'].between(inicio,fim)]
        if apenas:
            linhas = linhas.loc[linhas['Consenso']]
        if aa:
            linhas = linhas.loc[linhas['Aminoacido'].isin(aa)]
        st.caption(f'{len(linhas)} resíduos exibidos. O maior score ABCpred é o maior score dos peptídeos que cobrem a posição; não é uma predição por aminoácido do ABCpred.')
        tabela(linhas)
        st.download_button('Baixar aminoácidos exibidos (CSV)',csv_bytes(linhas),f'{proteina}_aminoacidos.csv','text/csv')
    with trechos:
        st.subheader('Trechos contínuos previstos pelas duas ferramentas')
        min_regiao = st.number_input('Comprimento mínimo da região (aa)',1,len(registro['sequencia']),1,key='reg_min_'+chave)
        regioes_visiveis = regioes.loc[regioes['Tamanho'] >= min_regiao]
        st.caption(f'{len(regioes_visiveis)} de {len(regioes)} regiões exibidas. Este filtro apenas oculta trechos; não une regiões separadas nem modifica a análise.')
        tabela(regioes_visiveis)
        st.download_button('Baixar regiões exibidas (CSV)',csv_bytes(regioes_visiveis),f'{proteina}_regioes.csv','text/csv')
    with origem:
        st.subheader('Parâmetros e dados utilizados')
        st.write('Os arquivos da execução selecionada permanecem intactos. Para obter outras janelas ou recuperar peptídeos abaixo do limiar ABCpred original, faça uma nova predição pelo main.py.')
        st.json({'origem':str(pasta),'predicao_original':{'abc':registro['abc'],'bepi':registro['bepi']},'analise_atual':resumo})
        st.download_button('Baixar análise atual (ZIP)',
                           pacote_analise(registro,pasta,residuos,peptideos,regioes,resumo,filtros,filtrados),
                           f'{proteina}_analise.zip','application/zip')
        st.download_button('Baixar gráfico atual (HTML)',figura.to_html(include_plotlyjs=True).encode('utf-8'),
                           f'{proteina}_grafico.html','text/html')
        with st.expander('Dados originais dos preditores'):
            st.write('ABCpred')
            tabela(abc)
            st.write('BepiPred')
            tabela(bepi)
        if problemas:
            with st.expander('Execuções ignoradas ou incompletas'):
                st.write(problemas)


if __name__ == '__main__':
    main()
