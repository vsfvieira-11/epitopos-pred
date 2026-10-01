import unittest
import pandas as pd
from comparacao import comparar


class ComparacaoTest(unittest.TestCase):
    def dados(self):
        seq = 'ACDEFGHIKL'
        abc = pd.DataFrame([[1, 'CDEF', 2, .8], [2, 'EFGH', 4, .7]],
                           columns=['Rank', 'Sequência', 'Posição de Início', 'Score'])
        bp = pd.DataFrame({'Posição de Início': range(1, 11), 'Aminoácido': list(seq),
                           'Score_Bruto_BP3': [.1, .8, .8, .1, .8, .8, .8, .1, .1, .8],
                           'Score_Linear_BP3': [.2]*10})
        return seq, abc, bp

    def test_overlap_counts_once_and_regions_have_inclusive_ends(self):
        seq, abc, bp = self.dados()
        residuos, peptideos, regioes, resumo = comparar(seq, abc, bp, .5)
        self.assertEqual(resumo['residuos_abc'], 6)
        self.assertEqual(resumo['residuos_consenso'], 5)
        self.assertEqual(regioes[['Inicio', 'Fim']].values.tolist(), [[2, 3], [5, 7]])
        self.assertEqual(regioes['Sequencia'].tolist(), ['CD', 'FGH'])
        self.assertEqual(peptideos['Residuos_BP_positivos'].tolist(), [3, 3])
        self.assertEqual(residuos.loc[4, 'ABC_n_peptideos'], 2)

    def test_linear_mode_is_separate(self):
        seq, abc, bp = self.dados()
        _, _, _, resumo = comparar(seq, abc, bp, .2, 'linear')
        self.assertEqual(resumo['residuos_bp'], 10)
        self.assertEqual(resumo['residuos_consenso'], 6)

    def test_empty_predictions_are_not_missing_results(self):
        seq, abc, bp = self.dados()
        _, peptideos, regioes, resumo = comparar(seq, abc.iloc[:0], bp, .9)
        self.assertTrue(peptideos.empty)
        self.assertTrue(regioes.empty)
        self.assertIsNone(resumo['jaccard_cobertura'])

    def test_wrong_sequence_and_positions_rejected(self):
        seq, abc, bp = self.dados()
        abc.loc[0, 'Sequência'] = 'AAAA'
        with self.assertRaises(ValueError):
            comparar(seq, abc, bp)
        seq, abc, bp = self.dados()
        bp.loc[0, 'Posição de Início'] = 2
        with self.assertRaises(ValueError):
            comparar(seq, abc, bp)
