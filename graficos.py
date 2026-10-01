"""Relatório Plotly independente, que pode ser aberto sem internet."""
from pathlib import Path
import html
import plotly.graph_objects as go
from plotly.subplots import make_subplots


def criar_figura(nome, residuos, resumo):
    pos = list(residuos['Posicao'])
    aa = list(residuos['Aminoacido'])
    fig = make_subplots(rows=3, cols=1, shared_xaxes=True,
                        row_heights=[.55, .25, .20], vertical_spacing=.09,
                        subplot_titles=['Scores BepiPred por aminoácido',
                                        'Cobertura e consenso (cor = presente)',
                                        'Quantidade de peptídeos ABCpred que cobrem cada posição'])
    for coluna, nome_score, cor, modo in [('Score_BP_bruto', 'Bruto', '#64748b', 'bruto'),
                                         ('Score_BP_linear', 'Linear suavizado', '#2563eb', 'linear')]:
        fig.add_trace(go.Scatter(x=pos, y=list(residuos[coluna]), mode='lines', name=nome_score,
                                line={'color': cor, 'width': 3 if resumo['score_bp_usado']==modo else 1.5},
                                customdata=aa, hovertemplate='Posição %{x} · %{customdata}<br>Score %{y:.4f}<extra>%{fullData.name}</extra>'), row=1, col=1)
    fig.add_hline(y=resumo['limiar_bp'], line_dash='dash', line_color='#dc2626',
                  annotation_text=f"Limiar {resumo['limiar_bp']:g} · classificação pelo score {resumo['score_bp_usado']}",
                  annotation_bgcolor='rgba(255,255,255,0.9)', annotation_font_size=10,
                  row=1, col=1)
    faixas = ['ABCpred', 'BepiPred ('+resumo['score_bp_usado']+')', 'Consenso']
    z = [[int(v) for v in residuos[c]] for c in ['ABC_coberto', 'BP_positivo', 'Consenso']]
    fig.add_trace(go.Heatmap(x=pos, y=faixas, z=z, zmin=0, zmax=1,
                            colorscale=[[0, '#f1f5f9'], [1, '#0f766e']], showscale=False,
                            customdata=[aa]*3,
                            hovertemplate='Posição %{x} · %{customdata}<br>%{y}: %{z} (0=ausente, 1=presente)<extra></extra>'), row=2, col=1)
    fig.add_trace(go.Bar(x=pos, y=list(residuos['ABC_n_peptideos']), name='Peptídeos ABCpred',
                        marker_color='#7c3aed', customdata=aa,
                        hovertemplate='Posição %{x} · %{customdata}<br>%{y} peptídeo(s)<extra></extra>'), row=3, col=1)
    fig.update_layout(template='plotly_white', height=900,
                      title={'text': f"{html.escape(nome)} · {resumo['residuos_consenso']}/{resumo['tamanho_proteina']} resíduos em consenso "
                                     f"({resumo['fracao_consenso']:.1%}) · {resumo['regioes_consenso']} regiões"},
                      margin={'l':140, 'r':40, 't':110, 'b':70},
                      legend={'orientation':'h', 'y':1.07, 'x':0}, hovermode='closest')
    fig.update_yaxes(title_text='Score BepiPred', range=[0, 1], row=1, col=1)
    fig.update_yaxes(autorange='reversed', row=2, col=1)
    fig.update_yaxes(title_text='Nº peptídeos', rangemode='tozero', row=3, col=1)
    fig.update_xaxes(range=[.5, len(pos)+.5])
    fig.update_xaxes(title_text='Posição na proteína (início em 1)', rangeslider_visible=True, row=3, col=1)
    return fig


def salvar_grafico(nome, residuos, resumo, destino):
    destino = Path(destino)
    fig = criar_figura(nome, residuos, resumo)
    fig.write_html(destino, include_plotlyjs=True, full_html=True,
                   config={'responsive': True, 'displaylogo': False,
                           'toImageButtonOptions': {'format':'svg', 'filename':nome+'_comparacao'}})
    return destino
