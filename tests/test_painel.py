import json
import os
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch
import pandas as pd
from painel_dados import listar_execucoes, carregar_execucao, analisar_local, filtrar_candidatos


class PainelTest(unittest.TestCase):
    def criar_dados(self, raiz):
        pasta = Path(raiz)/'Teste_execucao'
        pasta.mkdir()
        seq = 'ACDEFGHIKLMNPQRSTVWY'
        registro = {'proteina':'Teste', 'sequencia':seq,
                    'abc':{'limiar':.51,'janela':16,'filtro':'off'},
                    'bepi':{'score_comparacao':'linear','limiar_comparacao':.1512},
                    'estado':{'abc':'concluido','bepi':'concluido'}}
        (pasta/'execucao.json').write_text(json.dumps(registro),encoding='utf-8')
        pd.DataFrame([[1,seq[:16],1,.8],[2,seq[4:],5,.6]],
                     columns=['Rank','Sequência','Posição de Início','Score']).to_csv(pasta/'abc_resultado.csv',index=False)
        pd.DataFrame({'Posição de Início':range(1,21),'Aminoácido':list(seq),
                      'Score_Bruto_BP3':[.1]*20,'Score_Linear_BP3':[.2]*20}).to_csv(pasta/'bepi_resultado.csv',index=False)
        return pasta

    def test_local_cut_recomputes_and_candidate_filter_does_not(self):
        with tempfile.TemporaryDirectory() as tmp:
            pasta = self.criar_dados(tmp)
            self.assertEqual(len(listar_execucoes(tmp)[0]),1)
            registro,abc,bp = carregar_execucao(pasta)
            r,p,g,s = analisar_local(registro,abc,bp,'linear',.1512,.7)
            self.assertEqual(s['residuos_consenso'],16)
            self.assertEqual(len(p),1)
            self.assertEqual(len(filtrar_candidatos(p,fracao_min=1,tamanho_min=17)),0)
            self.assertEqual(s['residuos_consenso'],16)
            with self.assertRaises(ValueError):
                analisar_local(registro,abc,bp,'linear',.1512,.4)

    def test_app_reacts_to_score_without_network(self):
        from streamlit.testing.v1 import AppTest
        with tempfile.TemporaryDirectory() as tmp:
            self.criar_dados(tmp)
            with patch.dict(os.environ, {'EPITOPOS_RESULTADOS':tmp}), patch('requests.post') as post, patch('requests.get') as get:
                app = AppTest.from_file(str(Path(__file__).resolve().parents[1]/'app.py')).run(timeout=30)
                self.assertEqual(len(app.exception),0)
                self.assertEqual(app.metric[2].value,'20')
                app.radio[0].set_value('bruto').run(timeout=30)
                self.assertEqual(len(app.exception),0)
                self.assertEqual(app.metric[2].value,'0')
                post.assert_not_called()
                get.assert_not_called()
