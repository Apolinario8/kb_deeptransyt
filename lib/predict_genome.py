import argparse
import os
import gc
import time
import torch
import numpy as np
import pandas as pd
from Bio import SeqIO
from propythia.protein.sequence import ReadSequence
from data.DNN_binary import DNN_binary
from data.DNN import DNN
import logging
import esm
import json


logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

model_info = {}

def load_sequences(file_path: str) -> pd.DataFrame:
    sequences = []
    file_extension = file_path.split('.')[-1].lower()
    
    if file_extension in ['fasta', 'faa']:
        for record in SeqIO.parse(file_path, "fasta"):
            sequences.append((record.id, str(record.seq)))
    elif file_extension == 'txt':
        with open(file_path, 'r') as file:
            id = None
            seq = []
            for line in file:
                line = line.strip()
                if line.startswith('>'):
                    if id is not None:
                        sequences.append((id, ''.join(seq)))
                    id = line[1:].strip()
                    seq = []
                else:
                    seq.append(line)
            if id is not None:
                sequences.append((id, ''.join(seq)))
    else:
        raise ValueError("Unsupported file format. Supported formats are: .fasta, .faa and .txt")
    
    df = pd.DataFrame(sequences, columns=['ID', 'Sequence'])
    return df

def preprocess_sequences(df: pd.DataFrame, max_length: int = 600) -> pd.DataFrame:
    processed_sequences = []

    for _, row in df.iterrows():
        id, seq = row['ID'], row['Sequence']
        if len(seq) < 50:
            continue  # remove fragments
        if len(seq) < max_length:
            seq = seq.ljust(max_length, '-')  # Pad sequences shorter than max_length with '-'
        elif len(seq) > max_length:
            seq = seq[:max_length]  # Truncate sequences longer than max_length
        processed_sequences.append((id, seq))

    df = pd.DataFrame(processed_sequences, columns=['ID', 'Sequence'])

    read_seqs = ReadSequence()
    res = read_seqs.par_preprocessing(dataset=df, col='Sequence', B='N', Z='Q', U='C', O='K', J='I', X='')   # remove ambiguous aa

    return res

def create_encodings(df: pd.DataFrame, input_filename: str, model_name: str = 'esm2_t33_650M_UR50D', batch_size: int = 8, output_dir: str = 'encodings_genome', preprocess: bool = True) -> tuple:
    global model_info

    if preprocess:
        df = preprocess_sequences(df)

    repr_layer_map = {
        "esm2_t33": 33,
        "esm2_t6": 6,
        "esm2_t12": 12,
        "esm2_t30": 30,
        "esm2_t48": 48,
        "esm2_t36": 36
    }

    repr_layer = next((layer for name, layer in repr_layer_map.items() if name in model_name), None)
    if repr_layer is None:
        raise ValueError("Invalid model name provided")

    model, alphabet = esm.pretrained.load_model_and_alphabet(model_name)
    batch_converter = alphabet.get_batch_converter()

    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    model = model.to(device)

    labels = df["ID"].tolist()
    sequences = df["Sequence"].tolist()
    size = len(df["Sequence"][0])

    # Create batches
    data_batches = [(labels[start:start+batch_size], sequences[start:start+batch_size]) for start in range(0, len(labels), batch_size)]

    sequence_representations = []
    start_time = time.time()

    for batch_labels, batch_sequences in data_batches:
        batch_data = [(label, sequence) for label, sequence in zip(batch_labels, batch_sequences)]
        batch_labels, batch_strs, batch_tokens = batch_converter(batch_data)
        batch_tokens = batch_tokens.to(device)
        batch_lens = (batch_tokens != alphabet.padding_idx).sum(1)

        with torch.no_grad():
            results = model(batch_tokens, repr_layers=[repr_layer], return_contacts=True)

        token_representations = results["representations"][repr_layer].cpu()

        for i, tokens_len in enumerate(batch_lens):
            sequence_representations.append(token_representations[i, 1:tokens_len - 1].mean(0).cpu())     # Generate per-sequence representations

    total_time = time.time() - start_time
    logging.info(f"Total processing time: {total_time} seconds.")

    os.makedirs(output_dir, exist_ok=True)
    input_file_base = os.path.splitext(os.path.basename(input_filename))[0]
    encodings_filename = f"{output_dir}/{input_file_base}_encodings.npy"

    np.save(encodings_filename, np.array(sequence_representations))  

    model_info = {
        'model_name': model_name,
        'encodings_filename': encodings_filename,
        'aa_length': size,
    }

    # Clean up to free memory
    del model, alphabet, batch_converter, token_representations, results
    gc.collect()
    torch.cuda.empty_cache()

    return np.array(sequence_representations), labels

def predict_binary(encodings: np.ndarray, accession: list):
    model_path = 'models/DNN_allclasses.ckpt'

    model = DNN_binary.load_from_checkpoint(model_path)
    device = torch.device('cpu')
    model = model.to(device)
    model.eval()

    tensor_encodings = torch.tensor(encodings, dtype=torch.float32)
    with torch.no_grad():
        predictions = torch.sigmoid(model(tensor_encodings)).numpy().flatten()

    df_binary_predictions = pd.DataFrame({'Accession': accession, "Binary_Predictions": predictions})

    binary_labels = (predictions > 0.5).astype(int)
    num_transporters = np.sum(binary_labels)
    logging.info(f"Number of predicted transporters in the genome: {num_transporters}")

    return df_binary_predictions, binary_labels

def predict_family(transporter_encodings: np.ndarray, transporter_accessions: list) -> pd.DataFrame:    
    family_models = [
        {'path': '../data/models/family_DNN_no9_10.ckpt', 'num_classes': 402, 'column_name': 'PredictedFamily_>10', 'label_map': '../data/mappings/mapping_family10.json'},
        {'path': '../data/models/family_DNN_no9_12.ckpt', 'num_classes': 328, 'column_name': 'PredictedFamily_>12', 'label_map': '../data/mmappings/mapping_family12.json'},
        {'path': '../data/models/family_DNN_no9_15.ckpt', 'num_classes': 279, 'column_name': 'PredictedFamily_>15', 'label_map': '../data/mmappings/mapping_family15.json'},
        {'path': '../data/models/family_DNN_no9_20.ckpt', 'num_classes': 196, 'column_name': 'PredictedFamily_>20', 'label_map': '../data/mmappings/mapping_family20.json'},
        {'path': '../data/models/family_DNN_no9_30.ckpt', 'num_classes': 109, 'column_name': 'PredictedFamily_>30', 'label_map': '../data/mmappings/mapping_family30.json'},
        {'path': '../data/models/family_DNN_no9_40.ckpt', 'num_classes': 75, 'column_name': 'PredictedFamily_>40', 'label_map': '../data/mmappings/mapping_family40.json'},
        {'path': '../data/models/family_DNN_no9_50.ckpt', 'num_classes': 51, 'column_name': 'PredictedFamily_>50', 'label_map': '../data/mmappings/mapping_family50.json'},
    ]

    df_family_predictions = pd.DataFrame({'Accession': transporter_accessions})

    for model_info in family_models:
        num_classes = model_info['num_classes']
        
        family_model = DNN.load_from_checkpoint(checkpoint_path=model_info['path'], num_classes=num_classes)
        device = torch.device('cpu')
        family_model = family_model.to(device)
        family_model.eval()

        logging.info(f"Making Family Predictions for model with {num_classes} classes")
        transporter_tensor = torch.tensor(transporter_encodings, dtype=torch.float32)
        with torch.no_grad():
            transporter_predictions = family_model(transporter_tensor)
        predicted_families = transporter_predictions.argmax(dim=1).numpy()

        with open(model_info['label_map'], 'r') as f:
            label_map = json.load(f)

        predicted_family_labels = [label_map[str(label)] for label in predicted_families]

        df_family_predictions[model_info['column_name']] = predicted_family_labels

    return df_family_predictions

def predict_subfamily(transporter_encodings: np.ndarray, transporter_accessions: list) -> pd.DataFrame:
    subfamily_models = [
        {'path': '../data/models/subfamily_DNN_no9_15.ckpt', 'num_classes': 267, 'column_name': 'PredictedsubFamily_>15', 'label_map': '../data/mappings/mapping_subfamily15.json'},
        {'path': '../data/models/subfamily_DNN_no9_20.ckpt', 'num_classes': 159, 'column_name': 'PredictedsubFamily_>20', 'label_map': '../data/mappings/mapping_subfamily20.json'},
        {'path': '../data/models/subfamily_DNN_no9_30.ckpt', 'num_classes': 75, 'column_name': 'PredictedsubFamily_>30', 'label_map': '../data/mappings/mapping_subfamily30.json'},
        {'path': '../data/models/subfamily_DNN_no9_50.ckpt', 'num_classes': 30, 'column_name': 'PredictedSubFamily_>50', 'label_map': '../data/mappings/mapping_subfamily50.json'}
    ]

    df_subfamily_predictions = pd.DataFrame({'Accession': transporter_accessions})

    for model_info in subfamily_models:
        num_classes = model_info['num_classes']
        
        subfamily_model = DNN.load_from_checkpoint(checkpoint_path=model_info['path'], num_classes=num_classes)
        device = torch.device('cpu')
        subfamily_model = subfamily_model.to(device)
        subfamily_model.eval()

        logging.info(f"Making Sub-Family Predictions for model with {num_classes} classes")
        transporter_tensor = torch.tensor(transporter_encodings, dtype=torch.float32)
        with torch.no_grad():
            transporter_predictions = subfamily_model(transporter_tensor)
        predicted_subfamilies = transporter_predictions.argmax(dim=1).numpy()

        with open(model_info['label_map'], 'r') as f:
            label_map = json.load(f)

        predicted_subfamily_labels = [label_map[str(label)] for label in predicted_subfamilies]

        df_subfamily_predictions[model_info['column_name']] = predicted_subfamily_labels

    return df_subfamily_predictions

def predict_metabolic_important(transporter_encodings: np.ndarray, transporter_accessions: list) -> pd.DataFrame:     
    
    model_path = '..data/models/newdataset_test_no9.ckpt'
    label_map_path = '../data/mappings/mapping_newdataset.json'
    
    new_model = DNN.load_from_checkpoint(checkpoint_path=model_path, num_classes=4)
    device = torch.device('cpu')
    new_model = new_model.to(device)
    new_model.eval()

    logging.info(f"Making predictions for most important labels")
    transporter_tensor = torch.tensor(transporter_encodings, dtype=torch.float32)
    with torch.no_grad():
        transporter_predictions = new_model(transporter_tensor)
    predicted_labels = transporter_predictions.argmax(dim=1).numpy()

    with open(label_map_path, 'r') as f:
        label_map = json.load(f)

    predicted_labels_names = [label_map[str(label)] for label in predicted_labels]

    df_new_predictions = pd.DataFrame({'Accession': transporter_accessions, 'Predicted_newdataset': predicted_labels_names})

    logging.info(f"Predictions saved to results/metabolic_models_predictions.csv")
    return df_new_predictions

def process_genome(input_file: str, output_file: str, preprocess: bool):
    df_sequences = load_sequences(input_file)
    encodings, labels = create_encodings(df_sequences, input_file, preprocess=preprocess)
    
    df_binary_predictions, binary_labels = predict_binary(encodings, labels)
    
    positive_indices = np.where(binary_labels == 1)[0]
    positive_encodings = encodings[positive_indices]
    positive_labels = [labels[i] for i in positive_indices]

    df_family_predictions = predict_family(positive_encodings, positive_labels)
    df_subfamily_predictions = predict_subfamily(positive_encodings, positive_labels)
    df_metabolic_predictions = predict_metabolic_important(positive_encodings, positive_labels)
    
    df_merged = df_binary_predictions.merge(df_family_predictions, on='Accession', how='left')
    df_merged = df_merged.merge(df_subfamily_predictions, on='Accession', how='left')
    df_merged = df_merged.merge(df_metabolic_predictions, on='Accession', how='left')
    
    df_merged.to_csv(output_file, index=False)
    logging.info(f"All predictions saved to {output_file}")