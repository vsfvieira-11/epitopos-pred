import contextlib
import io
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch
import pandas as pd
import main


class MainTest(unittest.TestCase):
    def test_one_service_failure_preserves_other_result(self):
        abc = pd.DataFrame([[1, 'ACDEFGHIKL', 1, .8]],
                           columns=['Rank', 'Sequência', 'Posição de Início', 'Score'])
        with tempfile.TemporaryDirectory() as tmp, \
             patch('sys.argv', ['main.py', '--nome', 'Teste', '--sequencia', 'ACDEFGHIKL', '--janela-abc', '10', '--saida', tmp]), \
             patch.object(main, 'buscar_epitepos_abcpred', return_value=abc), \
             patch.object(main, 'buscar_epitopos_bepipred', side_effect=RuntimeError('Servidor indisponível')), \
             contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(io.StringIO()):
            self.assertEqual(main.main(), 1)
            pasta = next(Path(tmp).iterdir())
            self.assertTrue((pasta/'abc_resultado.csv').exists())
            self.assertFalse((pasta/'comparacao_residuos.csv').exists())
            registro = json.loads((pasta/'execucao.json').read_text(encoding='utf-8'))
            self.assertEqual(registro['estado']['comparacao'], 'nao_realizada')
            self.assertIn('bepi', registro['erros'])
