import unittest
from unittest.mock import patch, MagicMock
import pandas as pd
import abcpred_scrapper as abc


class ABCpredTest(unittest.TestCase):
    def consultar(self, **opcoes):
        resposta = MagicMock(text='<html></html>')
        tabela = pd.DataFrame([[1, 'ACDEFGHIKLMNPQRS', 1, .8]],
                              columns=['Rank', 'Sequence', 'Start position', 'Score'])
        with patch.object(abc.requests, 'post', return_value=resposta) as post, \
             patch.object(abc.pd, 'read_html', return_value=[tabela]):
            resultado = abc.buscar_epitepos_abcpred('Teste', 'ACDEFGHIKLMNPQRSTVWY', **opcoes)
        self.assertEqual(len(resultado), 1)
        return post.call_args.kwargs

    def test_default_payload(self):
        args = self.consultar()
        self.assertEqual(args['data']['Threshold'], '0.51')
        self.assertEqual(args['data']['window'], '16')
        self.assertEqual(args['data']['filter'], 'off')

    def test_custom_payload(self):
        args = self.consultar(limiar=.8, janela=20, filtro='on', timeout=30)
        self.assertEqual(args['data']['Threshold'], '0.8')
        self.assertEqual(args['data']['window'], '20')
        self.assertEqual(args['data']['filter'], 'on')
        self.assertEqual(args['timeout'], 30)

    def test_invalid_parameters_do_not_submit(self):
        with patch.object(abc.requests, 'post') as post:
            for opcoes in [{'limiar':float('nan')}, {'janela':15}, {'filtro':'talvez'}, {'timeout':-1}]:
                with self.assertRaises(ValueError):
                    abc.buscar_epitepos_abcpred('Teste', 'ACDEFGHIKLMNPQRSTVWY', **opcoes)
            post.assert_not_called()
