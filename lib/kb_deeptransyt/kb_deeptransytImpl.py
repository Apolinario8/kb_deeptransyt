# -*- coding: utf-8 -*-
#BEGIN_HEADER
import os
import yaml
import logging

from DeepTranSyT.sequence_processing import load_sequences, preprocess_sequences, create_encodings
from DeepTranSyT.make_predictions import predict_binary, predict_family, predict_subfamily, predict_metabolic_important
from DeepTranSyT import main as run_deeptransyt

from installed_clients.DataFileUtilClient import DataFileUtil
from installed_clients.KBaseDataObjectToFileUtilsClient import KBaseDataObjectToFileUtils
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

    def __init__(self, config):
        #BEGIN_CONSTRUCTOR
        self.workspaceURL = config['workspace-url']
        self.scratch = os.path.abspath(config['scratch'])
        self.callback_url = os.environ['SDK_CALLBACK_URL']
        logging.basicConfig(format='%(created)s %(levelname)s: %(message)s',
                            level=logging.INFO)
        #END_CONSTRUCTOR
        self.dfu = DataFileUtil(self.callback_url)
        pass

    def run_kb_deeptransyt(self, ctx, params):
        """
        This function runs the deeptransyt model, processes input sequences, and returns a KBaseReport.
        :param params: instance of mapping from String to unspecified object
        :returns: instance of type "ReportResults" -> structure: parameter
           "report_name" of String, parameter "report_ref" of String
        """

        if not isinstance(params['input_genome'], str) or not len(params['input_genome']):
            raise ValueError('Pass in a valid genome reference string')

        # setup
        with open("/kb/module/kbase.yml", 'r') as stream:
            data_loaded = yaml.safe_load(stream)
        version = str(data_loaded['module-version'])
        input_genome = params['input_genome']

        # create Util objects
        wsClient = workspaceService(self.workspaceURL, token=ctx['token'])
        object_to_file_utils = KBaseDataObjectToFileUtils(self.callback_url, token=ctx['token'])

        # get genomes
        genome_dir = os.path.join(self.scratch, 'genomes')
        os.mkdir(genome_dir)
        genome_info = wsClient.get_object_info_new({'objects': [{'ref': input_genome}]})[0]
        genome_input_type = genome_info[2]
        faa_locs = list()
        genome_ref_dict = {}
        if 'GenomeSet' in genome_input_type:
            genomeSet_object = wsClient.get_objects2({'objects': [{'ref': input_genome}]})['data'][0]['data']
            for ref_dict in genomeSet_object['elements'].values():
                genome_ref = ref_dict['ref']
                name = wsClient.get_object_info_new({'objects': [{'ref': genome_ref}]})[0][1]
                genome_ref_dict[name] = genome_ref
        else:
            genome_ref_dict[genome_info[1]] = input_genome
        for genome_name, genome_ref in genome_ref_dict.items():
            # this makes the names match if you are doing a genome or genomeSet
            faa_file = '%s.faa' % genome_name
            faa_object = object_to_file_utils.GenomeToFASTA({
                "genome_ref": genome_ref,
                "file": faa_file,
                "dir": genome_dir,
                "console": [],
                "invalid_msgs": [],
                'residue_type': 'protein',
                'feature_type': 'CDS',
                'record_id_pattern': '%%feature_id%%',
                'record_desc_pattern': '[%%genome_id%%]',
                'case': 'upper',
                'linewrap': 50
            })
            faa_locs.append(faa_object['fasta_file_path'])

        output_dir = os.path.join(self.scratch, 'results')
        os.makedirs(output_dir, exist_ok=True)
        
        # Run DeepTranSyT package
        preprocess = params.get('preprocess', True)
        gpu = params.get('gpu', 3)
        
        run_deeptransyt(input_file=genome_ref, output_dir=output_dir, preprocess=preprocess, gpu=gpu)

        # ctx is the context object
        # return variables are: output
        #BEGIN run_kb_deeptransyt
        report = KBaseReport(self.callback_url)
        report_info = report.create({'report': {'objects_created':[],
                                                'text_message': params['input_genome']},
                                                'workspace_name': params['workspace_name']})

        output = {
            'report_name': report_info['name'],
            'report_ref': report_info['ref'],
        }
        #END run_kb_deeptransyt

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
