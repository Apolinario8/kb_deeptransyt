# -*- coding: utf-8 -*-
#BEGIN_HEADER
import os
import shutil
import uuid

from predict_genome import *

from installed_clients.DataFileUtilClient import DataFileUtil
from installed_clients.KBaseReportClient import KBaseReport
from installed_clients.WorkspaceClient import Workspace as workspaceService
#END_HEADER


class kb_deeptransyt:
    '''
    Module Name:
    kb_deeptransyt

    Module Description:
    A KBase module: kb_deeptransyt
    '''

    ######## WARNING FOR GEVENT USERS ####### noqa
    # Since asynchronous IO can lead to methods - even the same method -
    # interrupting each other, you must be *very* careful when using global
    # state. A method could easily clobber the state set by another while
    # the latter method is running.
    ######################################### noqa
    VERSION = "0.0.1"
    GIT_URL = "git@github.com:Fxe/kb_deeptransyt.git"
    GIT_COMMIT_HASH = "5855499b63c264782e53366e8822e14f624d0992"

    #BEGIN_CLASS_HEADER
    #END_CLASS_HEADER

    # config contains contents of config file in a hash or None if it couldn't
    # be found
    def _get_input_file_ref_from_params(self, params):
        if 'input_genome' in params:
            return params['input_genome']
        else:
            if 'input_ws' not in params and 'input_file' not in params:
                raise ValueError('Either the "input_genome" field or the "input_ws" with "input_file" fields must be set.')
            return str(params['input_ws']) + '/' + str(params['input_file'])

    def create_report(self, token, ws, uuid_string, output_dir):
        output_files = [{'shock_id': self.dfu.file_to_shock({
                                'file_path': os.path.join(output_dir, 'final_predictions.csv'),
                                'make_handle': 0,
                                'pack': 'zip'})['shock_id'],
                        'name': 'final_predictions.csv',
                        'label': 'Predicted Results',
                        'description': 'CSV file containing the predicted results'}]

        report_params = {
            'file_links': output_files,
            'workspace_name': ws,
            'report_object_name': 'kb_deeptransyt_report_' + uuid_string
        }

        kbase_report_client = KBaseReport(self.callback_url, token=token)
        output = kbase_report_client.create_extended_report(report_params)
        return output

    def __init__(self, config):
        self.workspaceURL = config['workspace-url']
        self.scratch = os.path.abspath(config['scratch'])
        self.callback_url = os.environ['SDK_CALLBACK_URL']
        self.dfu = DataFileUtil(self.callback_url)
        pass

    def run_kb_deeptransyt(self, ctx, params):
        """
        This function runs the deeptransyt model, processes input sequences, and returns a KBaseReport.
        :param params: instance of mapping from String to unspecified object
        :returns: instance of type "ReportResults" -> structure: parameter
           "report_name" of String, parameter "report_ref" of String
        """

        token = ctx['token']
        uuid_string = str(uuid.uuid4())
        output_dir = os.path.join(self.scratch, uuid_string)
        os.mkdir(output_dir)

        input_dir = params['input_dir']
        preprocess = params.get('preprocess', True)

        # Preprocessing
        df_sequences = load_sequences(input_dir)
        encodings, labels = create_encodings(df_sequences, input_dir, preprocess=preprocess)
        
        # Binary prediction
        df_binary_predictions, binary_labels = predict_binary(encodings, labels)

        # Family prediction
        positive_indices = np.where(binary_labels == 1)[0]
        positive_encodings = encodings[positive_indices]
        positive_labels = [labels[i] for i in positive_indices]
        
        df_family_predictions = predict_family(positive_encodings, positive_labels)
        df_subfamily_predictions = predict_subfamily(positive_encodings, positive_labels)
        df_metabolic_predictions = predict_metabolic_important(positive_encodings, positive_labels)
        
        # Merging results
        df_merged = df_binary_predictions.merge(df_family_predictions, on='Accession', how='left')
        df_merged = df_merged.merge(df_subfamily_predictions, on='Accession', how='left')
        df_merged = df_merged.merge(df_metabolic_predictions, on='Accession', how='left')
        
        # Saving final csv
        final_output_path = os.path.join(output_dir, 'final_predictions.csv')
        df_merged.to_csv(final_output_path, index=False)

        output = self.create_report(token, params['workspace_name'], uuid_string, output_dir)

        shutil.rmtree(output_dir, ignore_errors=True)

        return [output]
    
    def status(self, ctx):
        #BEGIN_STATUS
        returnVal = {'state': "OK",
                     'message': "",
                     'version': self.VERSION,
                     'git_url': self.GIT_URL,
                     'git_commit_hash': self.GIT_COMMIT_HASH}
        #END_STATUS
        return [returnVal]
