import unittest
from io import BytesIO
from zipfile import ZipFile
from unittest.mock import MagicMock, patch

import requests
import bepipred as bp

SEQ = 'ACDEFGHIKL'


def response(text='', content=b'', status=200, state=None):
    r = requests.Response()
    r.status_code = status
    r.url = bp.CGI + '?jobid=ABC-12'
    r._content = content or text.encode()
    if state is not None:
        if state == 'finished':
            r._content = b'<a href="/results/bepipred3_results.zip">Download</a>'
        else:
            r._content = ("launchcheck('" + state + "','ABC-12')").encode()
    return r


def zip_result(seq=SEQ):
    buf = BytesIO()
    with ZipFile(buf, 'w') as z:
        z.writestr('results/raw_output.csv',
                   'Accession,Residue,BepiPred-3.0 score,BepiPred-3.0 linear epitope score\n' +
                   '\n'.join(f'protein,{aa},0.2,0.3' for aa in seq))
    return buf.getvalue()


class BepiPredTest(unittest.TestCase):
    def session(self):
        client = MagicMock()
        client.__enter__.return_value = client
        client.post.return_value = response()
        return client

    def test_multipart_and_download(self):
        client = self.session()
        client.get.side_effect = [response(state='finished'), response(content=zip_result())]
        with patch.object(bp.requests, 'Session', return_value=client):
            df = bp.buscar_epitopos_bepipred('protein', SEQ)
        self.assertEqual(list(df['Posição de Início']), list(range(1, 11)))
        self.assertEqual(df.attrs['job_id'], 'ABC-12')
        self.assertIn('files', client.post.call_args.kwargs)
        self.assertEqual(client.post.call_args.kwargs['data']['roll_mean'], 'no')
        self.assertIn('\r\n', client.post.call_args.kwargs['data']['fasta'])

    def test_connection_interruption_resumes_without_resubmitting(self):
        client = self.session()
        client.get.side_effect = [requests.ConnectionError('Disconnected'),
                                 response(state='finished'), response(content=zip_result())]
        with patch.object(bp.requests, 'Session', return_value=client), patch.object(bp.time, 'sleep'):
            bp.buscar_epitopos_bepipred('protein', SEQ)
        client.post.assert_called_once()

    def test_finished_without_output_fails_immediately(self):
        client = self.session()
        client.get.side_effect = [response(text='<p>ERROR: Could not find output file</p>')]
        with patch.object(bp.requests, 'Session', return_value=client):
            with self.assertRaisesRegex(RuntimeError, 'não disponibilizou'):
                bp.buscar_epitopos_bepipred('protein', SEQ)

    def test_html_200_is_not_zip(self):
        client = self.session()
        client.get.side_effect = [response(state='queued'), response(text='<html>Waiting</html>'),
                                 response(state='finished'), response(content=zip_result())]
        with patch.object(bp.requests, 'Session', return_value=client), patch.object(bp.time, 'sleep'):
            df = bp.buscar_epitopos_bepipred('protein', SEQ)
        self.assertEqual(len(df), 10)

    def test_resume_does_not_post(self):
        client = self.session()
        client.get.side_effect = [response(state='finished'), response(content=zip_result())]
        with patch.object(bp.requests, 'Session', return_value=client):
            bp.buscar_epitopos_bepipred('protein', SEQ, job_id='ABC-12')
        client.post.assert_not_called()

    def test_truncated_result_rejected(self):
        with self.assertRaisesRegex(ValueError, 'incompleto'):
            bp.ler_resultado(zip_result('AC'), SEQ)

    def test_fasta_and_validation(self):
        self.assertEqual(bp.preparar_sequencia('>protein\nACDEF\nGHIKL\n'), SEQ)
        for seq in ['>one\n'+SEQ+'\n>two\n'+SEQ, 'ACD1', 'A'*1024]:
            with self.assertRaises(ValueError):
                bp.preparar_sequencia(seq)
