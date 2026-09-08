from dataclasses import dataclass, field, fields, asdict
import difflib
import re
import os
from pathlib import Path
import datetime
from .exceptions import ProcessingError
import warnings
import logging
import json

_DATE_PATTERN = re.compile(r'(\d{4}-\d{2}-\d{2})')


# Configure logger with nice formatting
def _setup_logger():
    """Set up a logger with nice formatted output."""
    logger = logging.getLogger('de_quack')

    if not logger.handlers:
        logger.setLevel(logging.DEBUG)
        
        handler = logging.StreamHandler()
        handler.setLevel(logging.DEBUG)

        formatter = logging.Formatter(
            '%(levelname)-8s | %(name)s | %(message)s',
            datefmt='%Y-%m-%d %H:%M:%S'
        )
        handler.setFormatter(formatter)
        logger.addHandler(handler)
    
    return logger


logger = _setup_logger()



def _warn(message: str, category=UserWarning, stacklevel=2):
    logger.warning(message)


def _info(message: str):
    logger.info(message)


def _debug(message: str):
    logger.debug(message)


def _error(message: str):
    """Log an error message."""
    logger.error(message)





CORE_QUERIES = {
    'get_experiment': '''
                        SELECT 
                            e.experiment_id,
                            e.model,
                            e.date,
                            e.file,
                            e.experiment_name,
                            e.contrast,
                            e.annotation_version,
                            e.normalization,
                            e.other_info AS extra_info,
                        FROM experimental_data e
                        
                        WHERE 
                            ($1 IS NULL OR e.experiment_id = $1) AND
                            ($2 IS NULL OR e.date = $2) AND
                            ($3 IS NULL OR e.model LIKE '%' || $3 || '%') AND
                            ($4 IS NULL OR e.file LIKE '%' || $4 || '%') AND
                            ($5 IS NULL OR e.experiment_name LIKE '%' || $5 || '%') AND
                            ($6 IS NULL OR e.contrast LIKE '%' || $6 || '%') AND
                            ($7 IS NULL OR e.annotation_version LIKE '%' || $7 || '%') AND
                            ($8 IS NULL OR e.normalization LIKE '%' || $8 || '%')
                     ''',
   
    
    'find_delete_experiment': '''
                        SELECT experiment_id FROM experimental_data
                        WHERE
                            ($1 IS NULL OR experiment_id = $1) AND
                            ($2 IS NULL OR date = $2) AND
                            ($3 IS NULL OR model LIKE '%' || $3 || '%') AND
                            ($4 IS NULL OR file LIKE '%' || $4 || '%') AND
                            ($5 IS NULL OR experiment_name LIKE '%' || $5 || '%') AND
                            ($6 IS NULL OR contrast LIKE '%' || $6 || '%') AND
                            ($7 IS NULL OR annotation_version LIKE '%' || $7 || '%') AND
                            ($8 IS NULL OR normalization LIKE '%' || $8 || '%')
                        ''',
    'delete_gene_results': ''' 
                        DELETE FROM gene_results
                        WHERE
                            experiment_id = ANY($1)
                        ''',
    'delete_experiment': '''
                        DELETE FROM experimental_data
                        WHERE
                            experiment_id = ANY($1) RETURNING *
                        ''',
                        
                       
   
    'insert_de_arrow': '''INSERT INTO gene_results 
                            (experiment_id, gene_symbol, ensembl_id, log2fc, logCPM, pvalue, padj, stat, other_info) 
                            SELECT $1 AS experiment_id, 
                                        COALESCE(g.symbol, sym.symbol, df.gene_symbol) AS gene_symbol,
                                        COALESCE(g.ensembl_id, sym.ensembl_id, df.ensembl_id) AS ensembl_id,
                                        log2fc, 
                                        logCPM, 
                                        pvalue, 
                                        padj, 
                                        stat, 
                                        other_info 
                                    FROM df 
                                    LEFT JOIN genes g ON
                                        g.ensembl_id = df.ensembl_id 
                                    LEFT JOIN genes sym ON
                                        sym.symbol = df.gene_symbol OR
                                        (sym.prev_symbol IS NOT NULL AND df.gene_symbol IS NOT NULL AND
                                        list_contains(sym.prev_symbol, df.gene_symbol)) OR
                                        (sym.alias_symbol IS NOT NULL AND df.gene_symbol IS NOT NULL AND
                                        list_contains(sym.alias_symbol, df.gene_symbol))
                                    WHERE
                                     $2 IS NULL OR g.species = $2 OR sym.species = $2
                                    ''',
    'insert_de_arrows': '''INSERT INTO gene_results
                            (experiment_id, gene_symbol, ensembl_id, log2fc, logCPM, pvalue, padj, stat, other_info) 
                                SELECT  ids.new_id AS experiment_id, 
                                        COALESCE(g.symbol, sym.symbol, df.gene_symbol) AS gene_symbol,
                                        COALESCE(g.ensembl_id, sym.ensembl_id, df.ensembl_id) AS ensembl_id,
                                        log2fc, 
                                        logCPM, 
                                        pvalue, 
                                        padj, 
                                        stat, 
                                        other_info 
                                    FROM df 
                                    LEFT JOIN genes g ON
                                        g.ensembl_id = df.ensembl_id 
                                    LEFT JOIN genes sym ON
                                        sym.symbol = df.gene_symbol OR
                                        (sym.prev_symbol IS NOT NULL AND df.gene_symbol IS NOT NULL AND
                                        list_contains(sym.prev_symbol, df.gene_symbol)) OR
                                        (sym.alias_symbol IS NOT NULL AND df.gene_symbol IS NOT NULL AND
                                        list_contains(sym.alias_symbol, df.gene_symbol))
                                    LEFT JOIN id_map ids ON
                                        ids.old_id = df.experiment_id
                                    WHERE
                                     $1 IS NULL OR g.species = $1 OR sym.species = $1
                                    '''
    }


gene_columns = {
        'gene_symbol': [
            'gene_symbol', 'symbol', 'hgnc_symbol', 'genesymbol',
        ],
        'gene_name':['gene', 'gene_name', 'genename', 
            'name'],
        'ensembl_id': [
            'ensembl_id', 'ensembl', 'gene_id', 'geneid', 'ensembl_gene_id', 
            'target_id', 'feature_id'
        ],
        'pvalue': [
        'pvalue', 'p_value', 'p.value', 'p-value', 'pval', 'p.val', 'p_val'
        ],
        'padj': [
        'padj', 'p.adj', 'p_adj', 'p.adjusted', 'p_adjusted', 
        'fdr', 'qval', 'qvalue', 'q-value', 'adj.p.val', 'adj.p.value'
        ],
        'stat': [
            'stat', 'f', 'lr', 't'
        ],
        'log2fc': [
            'log2foldchange', 'log2fc', 'log2_fc', 'log2.fc', 'logfc', 'log_fc', 'log.fc'
            ],
        'base_mean': [
        'basemean', 'base_mean',  'aveexpr', 'tpm', 'fpkm'
        ],
        'logCPM':['logcpm', 'log_cpm', 'log.cpm']
        }


def _find_date_and_file(info):
    if isinstance(info, str):
        info=os.path.abspath(info)
    else:
        return datetime.datetime.now().strftime('%Y-%m-%d'), 'in-memory object'

    if os.path.isfile(info):
        match = _DATE_PATTERN.search(info) # first tries to find the date using regex on the file name
        file_path=os.path.abspath(info)
        if match:
            return datetime.datetime.strptime(match.group(1), '%Y-%m-%d'), file_path
        else: # if the search doesn't work, then uses the path library to find the time, or else uses none
            path=Path(file_path)
            if path.stat().st_mtime:
                return datetime.datetime.fromtimestamp(path.stat().st_mtime).strftime('%Y-%m-%d'), file_path
            else:
                return datetime.datetime.now().strftime('%Y-%m-%d'), file_path

def _process_metadata(metadata: object):
    if isinstance(metadata, dict) or isinstance(metadata, list):
        return metadata
    
    import polars as pl
    if isinstance(metadata, pl.DataFrame):
        return metadata.to_dicts()
    try:
        import pandas as pd
        if isinstance(metadata, pd.DataFrame):
            return metadata.to_dict(orient='records')
    except ImportError:
        pass
    
    if isinstance(metadata, str):
        if not os.path.isfile(metadata):
            raise FileNotFoundError(f'File not found: {metadata}')
        if metadata.lower().endswith('.json'):
            with open(metadata, 'r') as f:
                return json.load(f)
        else:
            with open(metadata, 'r') as f:
                line = f.readline()
                separator = '\t' if '\t' in line else (';' if ';' in line else ',')
            df = pl.read_csv(metadata, separator=separator, null_values = ['', 'NA', 'NaN', 'nan'])
            return df.to_dicts()

def _try_process_metadata(metadata: object):
    try:
        return _process_metadata(metadata)
    except (FileNotFoundError, json.JSONDecodeError, pl.exceptions.ComputeError) as e:
        raise ProcessingError(f"Error processing metadata: {e}")

@dataclass
class ExperimentMetadata:
    experiment_name: str = None
    date: str = None
    contrast: str = None
    file: str = None
    normalization: str = None
    model: str = None
    annotation_version: str = None
    
    additional_info: dict = field(default_factory=dict)



    def to_dict(self, metadata: dict, info):
        """Convert the dataclass to a dictionary, including additional_info."""
        field_names = [f.name for f in fields(self)]
        field_name_set = set(field_names)
        for key, value in metadata.items():
            cleaned_key = key.lower().strip().replace(' ', '_')
            if cleaned_key in field_name_set:
                setattr(self, cleaned_key, value)
            else:
                matches = difflib.get_close_matches(cleaned_key, field_names, n=1, cutoff=0.8)
                if matches:
                    setattr(self, matches[0], value)
                else:
                    self.additional_info[key] = value

        core_dict = {
            'experiment_name': self.experiment_name,
            'date': self.date,
            'contrast': self.contrast,
            'file': self.file,
            'normalization': self.normalization,
            'model': self.model,
            'annotation_version': self.annotation_version,
            'additional_info': dict(self.additional_info),
        }
        if core_dict['date'] is None:
            self.date, _ = _find_date_and_file(info)
            core_dict['date'] = self.date
        else:
            try:
                core_dict['date'] = datetime.datetime.strptime(core_dict['date'], '%Y-%m-%d').strftime('%Y-%m-%d')
            except ValueError:
                _warn(f"Provided date '{core_dict['date']}' is not in the correct format (YYYY-MM-DD). Defaulting to current date.")
                self.date, _ = _find_date_and_file(info)
                core_dict['date'] = self.date
        if core_dict['file'] is None:
            _, file_path = _find_date_and_file(info)
            core_dict['file'] = file_path
        if core_dict['experiment_name'] is None and 'file' in self.additional_info:
            self.experiment_name = core_dict['file']
            core_dict['experiment_name'] = self.experiment_name 
        for field_name in field_names:
            if core_dict[field_name] is None:
                _warn(f'{field_name} is not provided in the metadata. {field_name} has been defaulted to unknown')
                core_dict[field_name] = 'unknown'

        other_info = core_dict.pop("additional_info")

        return core_dict, json.dumps(other_info)


       

