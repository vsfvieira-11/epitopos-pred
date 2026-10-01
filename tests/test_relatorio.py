import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch
import pandas as pd
from main import recomparar
from graficos import criar_figura
from comparacao import comparar


class RelatorioTest(unittest.TestCase):
    def dados(self):
        seq = 'ACDEFGHIKLMNPQRS'
        abc = pd.DataFrame([[1, seq, 1, .8]], columns=['Rank','Sequência','Posição de Início','Score'])
        bp = pd.DataFrame({'Posição de Início':range(1,17), 'Aminoácido':list(seq),
                           'Score_Bruto_BP3':[.1]*16, 'Score_Linear_BP3':[.2]*16})
        return seq, abc, bp

    def test_local_mode_uses_saved_parameters_and_never_calls_services(self):
        seq, abc, bp = self.dados()
        with tempfile.TemporaryDirectory() as tmp:
            origem = Path(tmp)/'original'
            origem.mkdir()
            fonte = {'proteina':'Teste', 'sequencia':seq, 'abc':{'limiar':.51,'janela':16,'filtro':'off'},
                     'bepi':{'score_comparacao':'linear','limiar_comparacao':.1512},
                     'estado':{'abc':'concluido','bepi':'concluido'}}
            (origem/'execucao.json').write_text(json.dumps(fonte), encoding='utf-8')
            abc.to_csv(origem/'abc_resultado.csv', index=False)
            bp.to_csv(origem/'bepi_resultado.csv', index=False)
            snapshot = {p.name:p.read_bytes() for p in origem.iterdir()}
            with patch('main.buscar_epitepos_abcpred') as a, patch('main.buscar_epitopos_bepipred') as b:
                destino = recomparar(origem, Path(tmp)/'saidas')
                a.assert_not_called()
                b.assert_not_called()
            r = json.loads((destino/'execucao.json').read_text(encoding='utf-8'))
            self.assertEqual(r['resumo']['residuos_consenso'],16)
            self.assertEqual(r['resumo']['score_bp_usado'],'linear')
            self.assertEqual(snapshot, {p.name:p.read_bytes() for p in origem.iterdir()})
            self.assertIn('plotly.js', (destino/'grafico.html').read_text(encoding='utf-8'))
            outro = recomparar(origem, Path(tmp)/'saidas', score='bruto')
            r = json.loads((outro/'execucao.json').read_text(encoding='utf-8'))
            self.assertEqual(r['resumo']['residuos_consenso'],0)

    def test_graph_shows_the_same_classification_as_csv(self):
        seq, abc, bp = self.dados()
        residuos, _, _, resumo = comparar(seq, abc, bp, .1512,'linear')
        fig = criar_figura('Teste', residuos, resumo)
        self.assertEqual(list(fig.data[0].x), list(range(1,17)))
        self.assertEqual(list(fig.data[2].z[2]), [1]*16)
        self.assertEqual(fig.layout.shapes[0].y0,.1512)
