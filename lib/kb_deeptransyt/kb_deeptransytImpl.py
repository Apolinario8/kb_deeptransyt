# -*- coding: utf-8 -*-
#BEGIN_HEADER
import logging
import os

from predict_genome import *

from installed_clients.KBaseReportClient import KBaseReport
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
    def __init__(self, config):
        #BEGIN_CONSTRUCTOR
        self.callback_url = os.environ['SDK_CALLBACK_URL']
        self.shared_folder = config['scratch']
        logging.basicConfig(format='%(created)s %(levelname)s: %(message)s',
                            level=logging.INFO)
        #END_CONSTRUCTOR
        pass


    def run_kb_deeptransyt(self, ctx, params):
        """
        This example function accepts any number of parameters and returns results in a KBaseReport
        :param params: instance of mapping from String to unspecified object
        :returns: instance of type "ReportResults" -> structure: parameter
           "report_name" of String, parameter "report_ref" of String
        """

        input_dir = params['input_dir']
        output_dir = params['output_dir']
        preprocess = params.get('preprocess', True)

        # Running the deeptransyt process
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
        
        # Merging 
        df_merged = df_binary_predictions.merge(df_family_predictions, on='Accession', how='left')
        df_merged = df_merged.merge(df_subfamily_predictions, on='Accession', how='left')
        df_merged = df_merged.merge(df_metabolic_predictions, on='Accession', how='left')
        
        # Saving final csv
        final_output_path = os.path.join(output_dir, 'final_predictions.csv')
        df_merged.to_csv(final_output_path, index=False)
        
        # ctx is the context object
        # return variables are: output
        #BEGIN run_kb_deeptransyt
        report = KBaseReport(self.callback_url)
        report_info = report.create({'report': {'objects_created':[],
                                                'text_message': params['parameter_1']},
                                                'workspace_name': params['workspace_name']})

        # deep_transyt.run(xxxxx)
        print('works!', params)

        output = {
            'report_name': report_info['name'],
            'report_ref': report_info['ref'],
        }
        #END run_kb_deeptransyt

        # At some point might do deeper type checking...
        if not isinstance(output, dict):
            raise ValueError('Method run_kb_deeptransyt return value ' +
                             'output is not type dict as required.')
        # return the results
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
