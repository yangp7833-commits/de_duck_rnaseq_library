import json
import os
import re
import time
import uuid
from importlib import resources
import importlib
import hashlib
import duckdb
import polars as pl
from typing import Sequence, TypeAlias
from .core import DeQuackling, experiment_columns
from .exceptions import DeQuackError, ProcessingError, DuplicateGeneTableError
from .utilities import ExperimentMetadata, gene_columns,  CORE_QUERIES, _setup_logger, _try_process_metadata
logger = _setup_logger()




_core_queries = CORE_QUERIES
_GENE_ALIAS_TO_COLUMN = {
    alias.lower(): canonical
    for canonical, aliases in gene_columns.items()
    for alias in aliases + [canonical]
}
_gene_columns = ['gene_symbol', 'ensembl_id', 'log2fc', 'logCPM', 'pvalue', 'padj', 'stat']

ExperimentId: TypeAlias = int 
ExperimentMetadataField: TypeAlias = str | int | float | bool | None | dict[str, object] | list[object]
ExperimentMetadataRecord: TypeAlias = dict[str, ExperimentMetadataField]
ExperimentMetadataMap: TypeAlias = dict[ExperimentId, ExperimentMetadataRecord]
MetadataInput: TypeAlias = dict|list[dict]|str|pl.DataFrame| None
MetadataTableSource: TypeAlias = pl.DataFrame | pl.LazyFrame | str | dict[str, object] | list[dict[str, object]]


def get_unique_conn() -> str:
    unique_id = uuid.uuid4().hex[:8]
    return f':memory:de_quack_{unique_id}'

def _clean_df(df: pl.DataFrame) -> pl.DataFrame:
    df = df.with_columns(pl.col('gene_symbol').fill_null(pl.lit('')).alias('gene_symbol'),
                        pl.col('ensembl_id').fill_null(pl.lit('')).alias('ensembl_id'),
                        pl.col('log2fc').fill_null(pl.lit(0.0)).alias('log2fc'),
                        pl.col('pvalue').fill_null(pl.lit(1.0)).alias('pvalue'),
                        pl.col('logCPM').fill_null(pl.lit(0.0)).alias('logCPM'),
                        pl.col('other_info').fill_null(pl.lit('{}')).alias('other_info'),
                        pl.col('padj').fill_null(pl.lit(1.0)).alias('padj'))
    return df



def _to_polars_table(table: object, ignore_errors: bool = False) -> pl.DataFrame:
    """
    Convert various table-like objects to a polars DataFrame.
    Accepts polars DataFrames, polars LazyFrames, pandas DataFrames, file paths (CSV, TSV, Parquet), dictionaries, and lists of dictionaries.
    Attempts imports as needed.
    """
    if isinstance(table, pl.DataFrame):
        return table
    if isinstance(table, pl.LazyFrame):
        return table.collect()

    try:
        import pandas as pd
    except ImportError:
        pd = None

    if pd is not None and isinstance(table, pd.DataFrame):
        return pl.from_pandas(table, nan_to_null = True)

    if isinstance(table, str):
        path = os.path.abspath(table)
        if not os.path.exists(path):
            raise FileNotFoundError(f'File not found: {path}')
        if path.lower().endswith('.parquet'):
            return pl.read_parquet(path)
        if path.lower().endswith('xls') or path.lower().endswith('xlsx'):
            try:
                import openpyxl
            except ImportError:
                raise ImportError('openpyxl is required to read Excel files. Please install it with `pip install openpyxl`.')
            return pl.read_excel(path, sheet_name=0)
        try:
            with open(path, 'r', encoding='utf-8') as f:
                first_line = f.readline()
        except Exception as e:
            raise ProcessingError(f"Error occurred while reading the file: {e}")
        separator = '\t' if '\t' in first_line else (';' if ';' in first_line else ',')
        return pl.read_csv(path, separator=separator, null_values = ['', 'NA', 'NaN', 'nan'], truncate_ragged_lines=ignore_errors, quote_char = None)

    if isinstance(table, dict):
        return pl.DataFrame(table)

    if isinstance(table, list):
        return pl.DataFrame(table)

    return pl.DataFrame(table)


def _check_ids(ids: Sequence[ExperimentId]) -> None:
    """ 
    Ensures that the provided list of IDs is valid. Raises an error if any ID occurs multiple times or if any ID is a string.
    """
    for item in ids:
        if ids.count(item) > 1:
            raise DeQuackError(f'id {item} occurs multiple times in object ids')
        if isinstance(item, str):
            raise TypeError('String found in id list')


def _get_experiment_attribute(
    metadata: ExperimentMetadataMap,
    field: str,
    ids: ExperimentId | Sequence[ExperimentId],
) -> dict[ExperimentId, ExperimentMetadataField | None]:
    """
    Used for getting experiment attributes from the metadata dictionary. Returns a dictionary of {id: attribute} for each id in the provided list.
    Used in __init__ and _from_arrow calls.
    """
    if isinstance(ids, list):
        return {item: metadata.get(item, {}).get(field) for item in ids}
    return {ids: metadata.get(ids, {}).get(field)}


def _get_gene_table(species: str) -> object:
    """
    Retrieves the path of the gene table for the specified species.
    """
    folder = resources.files('de_quack') / 'gene_tables'
    return folder / f'{species}_genes.parquet'


def _heal_genes(df: pl.DataFrame, species: str) -> pl.DataFrame:
    """
    Heal gene symbols and Ensembl IDs in the provided DataFrame using the specefied gene table.
    Done with polars joining. If experiment id is present (as in DeArrows initiation), the column will be included in the select statement.
    """
    exp_id = 'experiment_id' in df.columns
    df = df.lazy()
    genes = pl.scan_parquet(_get_gene_table(species)).select('symbol', 'ensembl_id')
    df = df.with_columns(pl.col('gene_symbol').fill_null(pl.lit('')).alias('gene_symbol'))
    df = df.join(genes, on='ensembl_id', how='left')
    df = df.with_columns(pl.col('ensembl_id').fill_null(pl.lit('')).alias('ensembl_id'))
    df = df.join(genes, left_on='gene_symbol', right_on='symbol', how='left', suffix='_on_symbol')
    if exp_id:
        df = df.select(
            pl.col('experiment_id'),
            pl.coalesce(['symbol', 'gene_symbol']).alias('gene_symbol'),
            pl.coalesce(['ensembl_id', 'ensembl_id_on_symbol']).alias('ensembl_id'),
            pl.all().exclude(['symbol', 'ensembl_id_on_symbol', 'gene_symbol', 'ensembl_id', 'experiment_id'])
        ).collect()
    else:
        df = df.select(
            pl.coalesce(['symbol', 'gene_symbol']).alias('gene_symbol'),
            pl.coalesce(['ensembl_id', 'ensembl_id_on_symbol']).alias('ensembl_id'),
            pl.all().exclude(['symbol', 'ensembl_id_on_symbol', 'gene_symbol', 'ensembl_id', 'experiment_id'])
        ).collect()
    unmatched_count = df.filter((pl.col('gene_symbol') == '') & (pl.col('ensembl_id') == '')).height
    if unmatched_count > 0:
        logger.warning(f'{unmatched_count} genes could not be matched to the {species} gene table.')
    return df


def _order_columns(df: object, columns: dict[str, str] | None = None) -> pl.DataFrame:
    """
    Orders the columns of the provided DataFrame according to the expected gene columns and any additional columns.
    If the provided DataFrame is missing any of the expected gene columns, those columns will be added with null values.
    If the provided DataFrame has additional columns, those columns will be added to the end of the DataFrame in the order they appear in the provided DataFrame.
    If BaseMean is detected, it will be converted to logCPM and added to the DataFrame. THe original basemean will be kept in the additional information.
    """
    if columns is None:
        columns = {}
    df = _to_polars_table(df)
    rename_map = {
        column: _GENE_ALIAS_TO_COLUMN[column.lower()]
        for column in df.columns
        if column.lower() in _GENE_ALIAS_TO_COLUMN and column != 'experiment_id'
    }
    for key, value in columns.items():
        if key in df.columns and key not in rename_map.keys():
            if value not in gene_columns:
                raise DeQuackError(f'Column {value} is not a valid gene column')
            rename_map[key] = value
        else:
            raise DeQuackError(f'Column {key} is not in the dataframe or is already mapped')
    if rename_map:
        df = df.rename(rename_map)

    if 'gene_name' in df.columns:
        sample = df.select('gene_name').to_series().head(1).to_list()
        if sample and sample[0] is not None and re.match(r'^ENS', str(sample[0])):
            if 'ensembl_id' not in df.columns:
                df = df.rename({'gene_name': 'ensembl_id'})
        elif 'gene_symbol' not in df.columns:
            df = df.rename({'gene_name': 'gene_symbol'})

    expressions = []
    if 'experiment_id' in df.columns:
        expressions.append(pl.col('experiment_id'))
    

    for column in _gene_columns:
        if column in df.columns:
            expressions.append(pl.col(column))
        else:
            if column == 'logCPM' and 'base_mean' in df.columns:
                expressions.append(pl.col('base_mean').add(1).log(base = 2).alias('logCPM'))
            else:
                expressions.append(pl.lit(None).alias(column))

    extra_columns = [column for column in df.columns if column not in _gene_columns and column != 'experiment_id']
    if extra_columns:
        if len(extra_columns) == 1 and extra_columns[0] == 'other_info':
            expressions.append(pl.col('other_info'))
        else:
            struct_columns = [column for column in extra_columns if column != 'other_info']
            if 'other_info' in df.columns:
                struct_columns.append('other_info')
            expressions.append(pl.struct(struct_columns).struct.json_encode().alias('other_info'))
    elif 'other_info' in df.columns:
        expressions.append(pl.col('other_info'))
    else:
        expressions.append(pl.lit(None).alias('other_info'))

    return df.select(expressions)






class DeArrow:

    def __init__(
        self,
        info: object,
        experiment_id: ExperimentId | None = None,
        metadata: MetadataInput | None = None,
        heal_genes: bool = False,
        species: str | None = None,
        columns: dict[str, str] | None = None,
        ignore_errors: bool = False,
        **fields: ExperimentMetadataField,
    ) -> None:
        """
        Initialize a DeArrow object with the provided data and metadata.
        Parameters
        ----------
        info : object
            A table-like object containing gene expression data. Can be a polars DataFrame, polars LazyFrame, pandas DataFrame, file path (CSV, TSV, Parquet), dictionary, or list of dictionaries.
            The object should have columns for gene identifiers and expression values.
        experiment_id : ExperimentId, optional
            An integer ID for the experiment. If not provided, defaults to 1.
        metadata : MetadataInput, optional
            A dictionary containing metadata for the experiment. If not provided, metadata can be provided as keyword arguments. Metadata can also be provided as a string path to a JSON or CSV file, or as a polars or pandas DataFrame. 
            If no metadata is provided, an error will be raised.
        heal_genes : bool, optional
        If True, attempts to heal gene identifiers. Defaults to False.
        species : str, optional
            The species of the gene data. Defaults to 'human'.
        columns : dict[str, str], optional
            A dictionary mapping column names to their expected types.
        **fields : ExperimentMetadataField
            Additional metadata fields as keyword arguments.
        ignore_errors : bool, optional
            Option to ignore errors when processing CSVs and TSVs. This is done by turning on the truncate_ragged_lines option in polars. Default is False.
            This option is only available in DeArrow/DeArrows initialization and is not available in the add_experiment method. This is because the add_experiment method is primarily used for adding DeArrow or DeArrows objects, which are already processed and do not require this option.
            If you want to ignore errors when adding a new experiment, please initialize the DeArrow or DeArrows object with ignore_errors set to True before adding it to the existing DeArrow object.
        Raises if no metadata is provided or if it is not a dictionary.
        """
        if columns is None:
            columns = {}
        if metadata is None:
            metadata = fields
        metadata = _try_process_metadata(metadata)
        if isinstance(metadata, list) and len(metadata) > 1:
            logger.warning('Multiple metadata records provided. Only the first will be used for this DeArrow object.')
        metadata = metadata[0] if isinstance(metadata, list) else metadata
        if len(metadata) == 0:
            raise ProcessingError('No metadata provided in de_arrow object')
        
        if experiment_id is None:
            experiment_id = 1
        
        if not isinstance(metadata, dict):
            raise ProcessingError('Metadata must be a dictionary')
        self._table, self.experiment_metadata = self.__class__._to_de_arrow(info, metadata, columns, experiment_id, heal_genes = heal_genes, species = species, ignore_errors = ignore_errors)
        self.name = _get_experiment_attribute(self.experiment_metadata, 'experiment_name', experiment_id)
        self.annotation_version = _get_experiment_attribute(self.experiment_metadata, 'annotation_version', experiment_id)
        self.contrast = _get_experiment_attribute(self.experiment_metadata, 'contrast', experiment_id)
        self.file = _get_experiment_attribute(self.experiment_metadata, 'file', experiment_id)
        self.date = _get_experiment_attribute(self.experiment_metadata, 'date', experiment_id)
        self.model = _get_experiment_attribute(self.experiment_metadata, 'model', experiment_id)
        self.id = experiment_id
    
    @classmethod
    def _to_de_arrow(
        self,
        info: object,
        metadata: ExperimentMetadataRecord,
        columns: dict[str, str] | None,
        experiment_id: ExperimentId | None = None,
        heal_genes: bool = False,
        species: str | None = None,
        ignore_errors: bool = False,
    ) -> tuple[pl.DataFrame, ExperimentMetadataMap]:
        """
        Convert a table-like object and metadata into a DeArrow-compatible polars DataFrame and metadata map.
        First converts metadata into an ExperimentMetadata object, then converts the table-like object into a polars DataFrame, orders the columns, heals gene identifiers if specified, and returns the cleaned DataFrame along with the metadata map.
        """
        try:
            df = _to_polars_table(info, ignore_errors=ignore_errors)
        except Exception as e:
            raise ProcessingError(f"Error occurred while converting info to polars table: {e}")
        metadata_fields, other_info=ExperimentMetadata().to_dict(metadata, info)
        for key, value in json.loads(other_info).items():
            metadata_fields[key]=value
        metadata_fields={experiment_id: metadata_fields}
        df = _order_columns(df, columns)
        if heal_genes == True:
            if species is None:
                species = 'human'
            df = _heal_genes(df, species)
        if 'experiment_id' in df.columns:
            df = df.drop('experiment_id')
        df = df.insert_column(0, pl.lit(experiment_id).alias('experiment_id'))
       
        return _clean_df(df), metadata_fields

   

    @classmethod
    def _from_arrow(
        cls,
        df: pl.DataFrame,
        metadata: ExperimentMetadataMap,
        experiment_id: ExperimentId | None,
    ) -> "DeArrow":
        """
        Convert a nanoarrow array back into a de_arrow object.
        
        Parameters
        ----------
        df : polars.DataFrame
            A polars DataFrame to convert
        metadata : dict
            A nested dictionary containing metadata {ID: {field: value}}
        
        Returns
        -------
        de_arrow
            A de_arrow object with the provided array and metadata
        """
        instance = cls.__new__(cls)
        instance._table = df
        instance.experiment_metadata = metadata
        instance.name = _get_experiment_attribute(metadata, 'experiment_name', experiment_id)
        instance.annotation_version = _get_experiment_attribute(metadata, 'annotation_version', experiment_id)
        instance.contrast = _get_experiment_attribute(metadata, 'contrast', experiment_id)
        instance.file = _get_experiment_attribute(metadata, 'file', experiment_id)
        instance.date = _get_experiment_attribute(metadata, 'date', experiment_id)
        instance.model = _get_experiment_attribute(metadata, 'model', experiment_id)
        instance.id = experiment_id
        return instance

        
    
    def __setattr__(self, name: str, value: object) -> None:
        if name in ('_table', 'experiment_metadata', 'name', 'annotation_version', 'contrast', 'file', 'date', 'id', 'model'):
            object.__setattr__(self, name, value)
        else:
            table = object.__getattribute__(self, '_table')
            setattr(table, name, value)
    
    def __getitem__(self, key: object) -> object:
        return self._table[key]
       

    def __len__(self) -> int:
        return len(self._table)
    
    def __repr__(self) -> str:
        return repr(self._table)
    
    def __str__(self) -> str:
        return(str(self._table))

    

    def __getattr__(self, name: str) -> object:
        """
        Override the default attribute access to delegate to the underlying polars DataFrame.
        Will intercept method calls and return a new DeArrow object if the result is a polars DataFrame. With columns preserved.
        If a polars dataframe is returned and the columns are different, an error will be raised to prevent metadata loss.
        If the result is not a polars DataFrame, it will be returned as is along with the experiment metadata.
        """
        table = object.__getattribute__(self, '_table')
        attr = getattr(table, name)
        if not callable(attr):
            return attr

        def wrapper(*args: object, **kwargs: object) -> object:
            result = attr(*args, **kwargs)
            if isinstance(result, pl.DataFrame):
                if result.columns == table.columns:
                    return DeArrows._finalize_table(self, result)
                else:
                    raise DeQuackError('DeArrows object columns are immutable.')
            return result, self.experiment_metadata

        return wrapper

    
    

    
    def get_significant_genes(self, log2fc: float = 1, padj: float = 0.05, logCPM: float = 1) -> "DeArrow":
        """
        Gets all significant genes from the DeArrow object based on the provided log2fc, padj, and logCPM thresholds.
        """
        df = self._table.filter(
            (pl.col('log2fc').abs() >= log2fc) &
            (pl.col('padj') <= padj) &
            (pl.col('logCPM') >= logCPM)
        )
        return self._from_arrow(df, self.experiment_metadata, self.id)

    def get_downregulated(self, log2fc: float = -1, padj: float = 0.05, logCPM: float = 1) -> "DeArrow":
        """
        Gets all downregulated genes from the DeArrow object based on the provided log2fc, padj, and logCPM thresholds.
        """
        df = self._table.filter(
            (pl.col('log2fc') <= log2fc) &
            (pl.col('padj') <= padj) &
            (pl.col('logCPM') >= logCPM)
        )
        return self._from_arrow(df, self.experiment_metadata, self.id)

    def get_upregulated(self, log2fc: float = 1, padj: float = 0.05, logCPM: float = 1) -> "DeArrow":
        """
        Gets all upregulated genes from the DeArrow object based on the provided log2fc, padj, and logCPM thresholds.
        """
        frame = self._table.filter(
            (pl.col('log2fc') >= log2fc) &
            (pl.col('padj') <= padj) &
            (pl.col('logCPM') >= logCPM)
        )
        return self._from_arrow(frame, self.experiment_metadata, self.id)

    def set_id(self, experiment_id: ExperimentId) -> "DeArrow":
        """
        Set the ID of the DeArrow object. Will change the dataframe as well as the metadata and attributes.
        """
        if isinstance(experiment_id, int):
            frame = self._table.with_columns(pl.lit(experiment_id).alias('experiment_id'))
            metadata = {experiment_id: self.experiment_metadata[self.id]}
            return self._from_arrow(frame, metadata, experiment_id)
        raise DeQuackError(f'set_ids expected an integer, not {type(experiment_id)}')


    def get_gene(
        self,
        gene_symbol: str | None = None,
        ensembl_id: str | None = None,
        experiment_id: ExperimentId | None = None,
    ) -> "DeArrow":
        """
        Get a specific gene from the DeArrow object by gene symbol or Ensembl ID. If an experiment ID is provided, the gene will be filtered by that experiment ID as well.
        """
        expression = pl.lit(True)
        if experiment_id is not None and 'experiment_id' in self.columns:
            expression = expression & (pl.col('experiment_id') == experiment_id)
        if gene_symbol is not None:
            expression = expression & (pl.col('gene_symbol') == gene_symbol)
        if ensembl_id is not None:
            expression = expression & (pl.col('ensembl_id') == ensembl_id)
        df = self._table.filter(expression)
        return self._from_arrow(df, self.experiment_metadata, self.id)
        

    def add_experiment(
        self,
        data: object,
        metadata: MetadataInput | None = None,
        experiment_id: ExperimentId | None = None,
        **fields: ExperimentMetadataField,
    ) -> "DeArrows":
        """
        Add a new experiment to the DeArrow object. The new experiment can be provided as a DeArrow, DeArrows, or a polars or pandas DataFrame or a file path. The metadata for the new experiment can be provided as a dictionary or as keyword arguments. If no metadata is provided, the metadata from the existing DeArrow object will be used. 
        The new experiment will be added to the existing DeArrow object and a new DeArrows object will be returned.
        Uses functions outside of the class.
        """

        ids = [self.id]
        metadata = _try_process_metadata(metadata)
        if metadata is None and fields:
            metadata = fields 
        if isinstance(data, DeArrow):
            df, meta, new_ids = _add_experiment_arrow(self._table, self.experiment_metadata, ids, data, experiment_id)
        elif isinstance(data, DeArrows):
            df, meta, new_ids = _add_experiment_arrows(self._table, self.experiment_metadata, ids, data, experiment_id)
        else:
            df, meta, new_ids = _add_experiment_data(self._table, self.experiment_metadata, ids, data, metadata, experiment_id)
        return DeArrows._from_arrow(df, meta, new_ids)
    

    
    def insert(self, file: str, initialize_gene_table: bool = False, species: str | None = None) -> None:
        if not isinstance(file, str):
            raise TypeError(f'file must be a string, not {type(file)}')
        if not os.path.exists(os.path.abspath(file)):
            raise FileNotFoundError(f'File not found: {file}')
        file = os.path.abspath(file)
        df = self._table.select(pl.all().exclude(['experiment_id']))
        df = df.with_columns(
        pl.col("gene_symbol").replace("", None),
        pl.col("ensembl_id").replace("", None)
        )        
        metadata_fields = self.experiment_metadata[self.id]
        other_info = {}
        for key in metadata_fields.keys():
            if key not in ['model', 'date', 'file', 'experiment_name', 'contrast', 'annotation_version', 'normalization']:
                other_info[key] = metadata_fields.pop(key)
        db = DeQuackling(file).connect()
        if initialize_gene_table == True:
            species = species or 'human'
            try:
                db.initialize_gene_table(species)
            except DuplicateGeneTableError:
                raise DuplicateGeneTableError(f'Gene table for {species} already exists in database. Please use a different species or set initialize_gene_table to False')
        data_signature = hashlib.sha256(str(df).encode()).hexdigest()
        result = db.conn.execute(
            "INSERT INTO experimental_data (model, date, file, experiment_name, contrast, annotation_version, normalization, other_info, data_signature, duckdb_version, de_quack_version) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?) RETURNING experiment_id",
            (metadata_fields.get('model'), metadata_fields.get('date'), metadata_fields.get('file'), metadata_fields.get('experiment_name'), metadata_fields.get('contrast'), metadata_fields.get('annotation_version'), metadata_fields.get('normalization'), other_info, data_signature, importlib.metadata.version('duckdb'), importlib.metadata.version('de_quack'))
        ).fetchall()
        db.conn.execute(_core_queries['insert_de_arrow'], (result[0][0], species))
        db.conn.commit()
        db.conn.close()
        

        
    def df(self) -> pl.DataFrame:
        return self._table.clone()
    
    def write_parquet(self, output_path: str, compression: str = 'zstd', compression_level: int = 3) -> None:
        """
        Calls the _write_parquet function to write the DeArrow object to parquet.
        """
        _write_parquet(self._table, self.experiment_metadata, output_path, compression=compression, compression_level=compression_level)


class DeArrows:
    def __new__(
        cls,
        *args: object,
        columns: dict[str, str] | None = None,
        metadata: list[ExperimentMetadataRecord] | ExperimentMetadataRecord | None = None,
        ids: list[ExperimentId] | None = None,
        keep_ids: bool = False,
        heal_genes: bool = False,
        species: str = 'human',
        ignore_errors: bool = False,
    ) -> object:
        """
        Checks if there are multiple metadata and provided tables. If there is only one metadata and one table, it will return a DeArrow object instead of a DeArrows object.
        """
        if columns is None:
            columns = {}
        if metadata is None:
            return super().__new__(cls)
        _try_process_metadata(metadata)
        if len(metadata) != len(args):
            if len(args) == 1 and isinstance(metadata, dict):
                return DeArrow(args[0], metadata=metadata)
        if len(args) == 1:
            return DeArrow(args[0], metadata=metadata[0])
        return super().__new__(cls)
    
    def __init__(
        self,
        *args: object,
        columns: dict[str, str] | None = None,
        metadata: list[ExperimentMetadataRecord] | ExperimentMetadataRecord | None = None,
        ids: list[ExperimentId] | None = None,
        heal_genes: bool = False,
        species: str = 'human',
        keep_ids: bool = False,
        ignore_errors: bool = False,
    ) -> None:
        if columns is None:
            columns = {}
        """Initialize a DeArrows object with multiple differential expression tables.
        
        Parameters
        ----------
        *args : str or file paths
            Variable length argument list of table references (file paths or identifiers)
        metadata : MetadataInput, optional
            Can be a list of dictionaries, a string path to a JSON or CSV file, or a polars or pandas DataFrame. Each dictionary should contain metadata for the corresponding table in args. If no metadata is provided, an error will be raised.
            Metadata for each table. Note that since this is DeArrows, the top JSON object in the JSON file must be a JSON array.
        ids : list of int, optional
            List of unique identifiers for each table. If not provided, will be generated automatically.
        heal_genes : bool, optional
            Whether to heal gene symbols and Ensembl IDs using the species-specific gene table. Default is
        columns: dict, optional
            A dictionary mapping column names in the input tables to the expected column names in the DeArrows object. Default is None.
            Use this if your input table has different columns than the implemented ones.
        keep_ids : bool, optional
            Whether to keep the ids in provided DeArrow objects. If False, new ids will be expected to be provided for DeArrow objects.
            If no ids are provided, new ids will be generated automatically. Default is False.
        species : str, optional
            The species to use for gene healing. Default is 'human'. If heal_genes is True, this parameter is required to specify the species for gene healing. If heal_genes is False, this parameter is ignored.
        ignore_errors : bool, optional
            Option to ignore errors when processing CSVs and TSVs. This is done by turning on the truncate_ragged_lines option in polars. Default is False.
        
        If a provided table already has an experiment_id column and it is NOT a DeArrow or DeArrows object, the id will be discarded. This is because maintaing equal lengths between the id and tables with keep_ids on is not possible when the id is already present in the table. 
        If you want to keep the id, please provide a DeArrow or DeArrows object instead of a table with an experiment_id column.
        """

        
        none_de_arrow = [table for table in args if not isinstance(table, DeArrow) and not isinstance(table, DeArrows)]
        if metadata is None:
            if len(none_de_arrow) == 0:
                metadata = []
            else:
                raise DeQuackError('Metadata must be provided for non-de_arrow tables')
        if len(args) < 1:
            raise DeQuackError('At least one table must be provided')
        if not isinstance(metadata, list):
            if len(none_de_arrow) == 1:
                metadata = [metadata]
            else:
                raise DeQuackError('Metadata must be provided as a list of dictionaries')
        if len(metadata) != len(none_de_arrow):
            raise DeQuackError('Number of metadata dictionaries does not match number of tables')
        if keep_ids == True:
            if ids is None:
                ids = [i + 1 for i in range(len(none_de_arrow))]
            if len(ids) != len(none_de_arrow):
                raise DeQuackError('Number of ids provided does not match number of tables')
            if isinstance(ids, int) and len(none_de_arrow) == 1:
                ids = [ids]
        else:
            if ids is None:
                if any(isinstance(arg, DeArrows) for arg in args):
                    raise DeQuackError('ids must be provided when combining DeArrows objects')
                ids = [i + 1 for i in range(len(args))]
            if len(ids) != len(args):
                raise DeQuackError('Number of ids provided does not match number of tables')
            
        _check_ids(ids)
        
        
        self._table, self.experiment_metadata, ids = self.__class__._from_tables(*args, metadata=metadata, ids=ids, columns=columns, heal_genes = heal_genes, species = species, keep_ids = keep_ids, ignore_errors = ignore_errors)
        self.id = ids
        self.name = _get_experiment_attribute(self.experiment_metadata, 'experiment_name', ids)
        self.annotation_version = _get_experiment_attribute(self.experiment_metadata, 'annotation_version', ids)
        self.contrast = _get_experiment_attribute(self.experiment_metadata, 'contrast', ids)
        self.file = _get_experiment_attribute(self.experiment_metadata, 'file', ids)
        self.date = _get_experiment_attribute(self.experiment_metadata, 'date', ids)
        self.model = _get_experiment_attribute(self.experiment_metadata, 'model', ids)

    
    def __getattr__(self, name: str) -> object:
        """
        Override the default attribute access to delegate to the underlying polars DataFrame.
        Will intercept method calls and return a new DeArrows object if the result is a polars DataFrame. With columns preserved.
        If a polars dataframe is returned but the columns are different, an error will be raised.
        This is done to ensure that metadata isn't accidentally dropped when performing operations on the underlying DataFrame.
        """
        table = object.__getattribute__(self, '_table')
        attr = getattr(table, name)
        if not callable(attr):
            return attr

        def wrapper(*args: object, **kwargs: object) -> object:
            result = attr(*args, **kwargs)
            if isinstance(result, pl.DataFrame):
                if result.columns == table.columns:
                    return self.__class__._finalize_table(self, result)
                else:
                    raise DeQuackError('DeArrows object columns are immutable.')
            return result, self.experiment_metadata

        return wrapper
    
    def __setattr__(self, name: str, value: object) -> None:
        if name in ('_table', 'experiment_metadata', 'name', 'annotation_version', 'contrast', 'file', 'date', 'id', 'model'):
            object.__setattr__(self, name, value)
        else:
            table = object.__getattribute__(self, '_table')
            setattr(table, name, value)
    
    def __getitem__(self, key: object) -> object:
        return self._table[key]
    
    def __len__(self) -> int:
        return len(self._table)
    
    def __repr__(self) -> str:
        return repr(self._table)
    
    def __str__(self) -> str:
        return(str(self._table))

    @classmethod
    def _from_tables(
        cls,
        *args: object,
        columns: dict[str, str] | None = {},
        metadata: list[ExperimentMetadataRecord] | None = None,
        ids: list[ExperimentId] | None = None,
        heal_genes: bool = False,
        species: str = 'human',
        keep_ids: bool = False,
        ignore_errors: bool = False,
    ) -> tuple[pl.DataFrame, ExperimentMetadataMap, list[ExperimentId]]:
        table = pl.DataFrame(schema = {'experiment_id': pl.Int32(), 'gene_symbol': pl.String(), 'ensembl_id': pl.String(), 'log2fc': pl.Float64(), 'logCPM': pl.Float64(), 'pvalue': pl.Float64(), 'padj': pl.Float64(), 'stat': pl.Float64(), 'other_info': pl.String()})
        frames = []
        meta_by_id = {}
        flat_ids = []
        meta_iter, ids_iter = iter(metadata or [{}]), iter(ids or [])
        for arg in args:
            if isinstance(arg, DeArrow):
                new_id = arg.id if keep_ids == True else next(ids_iter)
                arg = arg.with_columns(pl.lit(new_id).alias('experiment_id'))
                frames.append(arg._table)
                meta_by_id[new_id] = arg.experiment_metadata[arg.id]
                flat_ids.append(new_id)
            elif isinstance(arg, DeArrows):
                new_ids = arg.id if keep_ids == True else next(ids_iter)
                if not isinstance(new_ids, list) and not isinstance(new_ids, tuple) or len(new_ids) != len(arg.id):
                    raise DeQuackError('ids must be a list or tuple for DeArrows objects and must be of same id length for DeArrows objects')
                if not keep_ids:
                    arg = arg.with_columns(pl.col('experiment_id').replace(dict(zip(arg.id, new_ids))).alias('experiment_id'))
                frames.append(arg._table)
                flat_ids.extend(new_ids)
                for old_id, new_id in zip(arg.id, new_ids):
                    meta_by_id[new_id] = arg.experiment_metadata[old_id]
            else:
                new_id = next(ids_iter)
                try:
                    arg = _to_polars_table(arg, ignore_errors=ignore_errors)
                except Exception as e:
                    raise ProcessingError(f"Error occurred while converting info to polars table: {e}")
                arg = _order_columns(arg, columns)
                if 'experiment_id' in arg.columns:
                    arg = arg.drop('experiment_id')
                arg = arg.insert_column(0, pl.lit(new_id).alias('experiment_id'))
                flat_ids.append(new_id)
                meta = next(meta_iter)
                meta, other_info = ExperimentMetadata().to_dict(meta, arg)
                meta.update(json.loads(other_info))
                meta_by_id[new_id] = meta
                try:
                    frames.append(arg)
                except pl.exceptions.SchemaError as e:
                    raise DeQuackError(f'A table either has invalid data types or is missing required columns: {e}')
        table = pl.concat(frames, how = 'diagonal')
        if heal_genes == True:
            table = _heal_genes(table, species)
        _check_ids(flat_ids)
        return _clean_df(table), meta_by_id, flat_ids

    @classmethod
    def _from_arrow(
        cls,
        table: pl.DataFrame,
        metadata: ExperimentMetadataMap,
        ids: list[ExperimentId],
    ) -> "DeArrows":

        """ 
        Wrap a polars Dataframe back into a DeArrows object with the provided metadata and ids.
        """
        instance = object.__new__(cls)
        instance._table = table
        instance.experiment_metadata = metadata
        instance.id = ids
        instance.name = _get_experiment_attribute(metadata, 'experiment_name', ids)
        instance.annotation_version = _get_experiment_attribute(metadata, 'annotation_version', ids)
        instance.contrast = _get_experiment_attribute(metadata, 'contrast', ids)
        instance.file = _get_experiment_attribute(metadata, 'file', ids)
        instance.date = _get_experiment_attribute(metadata, 'date', ids)
        instance.model = _get_experiment_attribute(metadata, 'model', ids)
        return instance

    def set_id(self, ids: dict[ExperimentId, ExperimentId] | list[ExperimentId]) -> "DeArrows":
        """Set the IDs of the experiments in the DeArrows object."""
        if isinstance(ids, dict):
            return self._set_id_dict(ids)
        elif isinstance(ids, list):
            if len(ids) != len(self.id):
                raise DeQuackError('Number of provided ids does not match number of existing ids')
            return self._set_id_list(ids)
        raise DeQuackError(f'set_ids expected a dictionary or list, not {type(ids)}')
    
    def _set_id_dict(self, ids: dict[ExperimentId, ExperimentId]) -> "DeArrows":
        df = self._table
        for old_id, new_id in ids.items():
            if old_id not in self.id:
                raise DeQuackError(f'{old_id} is not in existing ids. Existing ids are {self.id}')
        _check_ids(list(ids.values()) + self.id)
        df = df.with_columns(pl.col('experiment_id').replace(ids).alias('experiment_id'))
        new_ids = self.id.copy()
        metadata = self.experiment_metadata.copy()
        for old_id, new_id in ids.items():
            index = self.id.index(old_id)
            new_ids[index] = new_id
            metadata[new_id] = self.experiment_metadata.get(old_id, {})

        return self._from_arrow(df, metadata, ids = new_ids)

    def _set_id_list(self, ids: list[ExperimentId]) -> "DeArrows":
        if len(ids) != len(self.id):
            raise DeQuackError('Number of provided ids does not match number of existing ids')
        _check_ids(ids)
        mapping = dict(zip(self.id, ids))
        df = self._table.with_columns(pl.col('experiment_id').replace(mapping).alias('experiment_id'))
        metadata = {new_id: self.experiment_metadata.get(old_id, {}) for old_id, new_id in mapping.items()}
        return self._from_arrow(df, metadata, ids=ids)
    
    def get_experiment(
        self,
        experiment_id: ExperimentId | list[ExperimentId] | None = None,
        name: str | list[str] | None = None,
        model: str | list[str] | None = None,
        annotation_version: str | list[str] | None = None,
        normalization: str | list[str] | None = None,
        date: str | None = None,
        contrast: str | list[str] | None = None,
        file: str | list[str] | None = None,
    ) -> "DeArrow | DeArrows":
        """
        Gets experiments from the DeArrows object based on the provided metadata. If no metadata is provided, all experiments will be returned. If multiple experiments are found, a DeArrows object will be returned. 
        If only one experiment is found, a DeArrow object will be returned.
        Uses the attributes from experiment metadata and gets the id to filter the dataframe. If no experiments are found, an error will be raised.
        """
        

        selected_ids = set(self.id)
        if experiment_id is not None:
            if isinstance(experiment_id, list):
                selected_ids = {item for item in self.id if item in experiment_id}
            else:
                selected_ids = {item for item in self.id if item == experiment_id}


            def matches(meta: ExperimentMetadataRecord) -> bool:
                def text_match(field_value: object, expected: object) -> bool:
                    if expected is None:
                        return True
                    if field_value is None:
                        return False
                    if isinstance(expected, list):
                        expected_values = [str(item).lower().strip() for item in expected]
                        return str(field_value).lower().strip() in expected_values
                    return str(expected).lower().strip() in str(field_value).lower().strip()

                return (
                    text_match(meta.get('experiment_name'), name) and
                    text_match(meta.get('model'), model) and
                    text_match(meta.get('annotation_version'), annotation_version) and
                    text_match(meta.get('normalization'), normalization) and
                    text_match(meta.get('contrast'), contrast) and
                    text_match(meta.get('file'), file) and
                    (date is None or str(meta.get('date')) == str(date))
                )

            selected_ids = {item for item in selected_ids if matches(self.experiment_metadata.get(item, {}))}

        if not selected_ids:
            raise DeQuackError('No experiments found with the provided metadata')

        frame = self._table.filter(pl.col('experiment_id').is_in(list(selected_ids))) if 'experiment_id' in self.columns else self._frame
        return self._finalize_table(frame)
    
    def _finalize_table(self, df: pl.DataFrame) -> "DeArrow | DeArrows":
        """
        Selects all unique experiment IDs from the provided DataFrame and retrieves the corresponding metadata. If no experiment IDs are found, a DeArrow object is returned with no ID. If one experiment ID is found, a DeArrow object is returned with that ID. 
        If multiple experiment IDs are found, a DeArrows object is returned with those IDs.
        Also changes the metadata and attributes of the returned object to match the experiment IDs.
        """
        df_ids = df.select('experiment_id').unique().to_series().to_list() if 'experiment_id' in df.columns else []
        metadata = {}
        for experiment_id in df_ids:
            metadata[experiment_id] = self.experiment_metadata.get(experiment_id, {})
        if len(df_ids) == 0:
            return DeArrow._from_arrow(df, metadata, None)
        if len(df_ids) == 1:
            return DeArrow._from_arrow(df, metadata, df_ids[0])
        return self._from_arrow(df, metadata, df_ids)
    
    def get_significant_genes(self, log2fc: float = 1, padj: float = 0.05, logCPM: float = 1) -> "DeArrows":
        df = self._table.filter(
            (pl.col('log2fc').abs() >= log2fc) &
            (pl.col('padj') <= padj) &
            (pl.col('logCPM') >= logCPM)
        )
        return self._finalize_table(df)

    def get_downregulated(self, log2fc: float = -1, padj: float = 0.05, logCPM: float = 1) -> "DeArrows":
        df = self._table.filter(
            (pl.col('log2fc') <= log2fc) &
            (pl.col('padj') <= padj) &
            (pl.col('logCPM') >= logCPM)
        )
        return self._finalize_table(df)

    def get_upregulated(self, log2fc: float = 1, padj: float = 0.05, logCPM: float = 1) -> "DeArrows":
        df = self._table.filter(
            (pl.col('log2fc') >= log2fc) &
            (pl.col('padj') <= padj) &
            (pl.col('logCPM') >= logCPM)
        )
        return self._finalize_table(df)
    
    def get_gene(
        self,
        gene_symbol: str | None = None,
        ensembl_id: str | None = None,
        experiment_id: ExperimentId | None = None,
    ) -> "DeArrows":
        expression = pl.lit(True)
        if experiment_id is not None and 'experiment_id' in self.columns:
            expression = expression & (pl.col('experiment_id') == experiment_id)
        if gene_symbol is not None:
            expression = expression & (pl.col('gene_symbol') == gene_symbol)
        if ensembl_id is not None:
            expression = expression & (pl.col('ensembl_id') == ensembl_id)
        df = self._table.filter(expression)
        return self._finalize_table(df)
        
    
    def insert(self, file: str, initialize_gene_table: bool = False, species: str | None = None) -> None:
        """
        Inserts the dataframe into a DuckDB file. If the file does not exist, it will be created. If the file exists, the dataframe will be appended to the existing data. If initialize_gene_table is True, a gene table will be created in the database for the specified species. 
        If the gene table already exists, an error will be raised. The species parameter is required if initialize_gene_table is True.
        Replaces empty strings in gene_symbol and ensembl_id columns with None to avoid issues with the database. Also calculates a data signature for the dataframe to ensure data integrity.
        Loops through all IDs and inserts them as experiments first. Then, retrieves new IDs generated from the database primary key sequence and creates an id_map for insertion of the table.
        A join will be performed on the id_map to replace old IDs with new IDs in the dataframe before insertion into the database. Finally, the dataframe will be inserted into the database and the connection will be closed.
        IDs in the DeArrows object will not be preserved, as there is a primary key constraint on all de_quack files.
        """
        if not isinstance(file, str):
            raise TypeError(f'file must be a string, not {type(file)}')
        if not os.path.exists(os.path.abspath(file)):
            raise FileNotFoundError(f'File not found: {file}')
        file = os.path.abspath(file)
        df = self._table
        df = df.with_columns(
        pl.col("gene_symbol").replace("", None),
        pl.col("ensembl_id").replace("", None)
        )
        db = DeQuackling(file).connect()
        old_ids = self.id
        new_ids = []
        for experiment_id in old_ids:        
            metadata_fields = self.experiment_metadata[experiment_id]
            other_info = {}
            for key in metadata_fields.keys():
                if key not in ['model', 'date', 'file', 'experiment_name', 'contrast', 'annotation_version', 'normalization']:
                    other_info[key] = metadata_fields.pop(key)
            data_signature = hashlib.sha256(str(df).encode()).hexdigest()
            result = db.conn.execute(
                "INSERT INTO experimental_data (model, date, file, experiment_name, contrast, annotation_version, normalization, other_info, data_signature, duckdb_version, de_quack_version) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?) RETURNING experiment_id",
                (metadata_fields.get('model'), metadata_fields.get('date'), metadata_fields.get('file'), metadata_fields.get('experiment_name'), metadata_fields.get('contrast'), metadata_fields.get('annotation_version'), metadata_fields.get('normalization'), other_info, data_signature, importlib.metadata.version('duckdb'), importlib.metadata.version('de_quack'))
            ).fetchall()
            new_ids.append(result[0][0])
        id_map = pl.DataFrame({'old_id': old_ids, 'new_id': new_ids}, schema = {'old_id': pl.Int32(), 'new_id': pl.Int32()})
        if initialize_gene_table == True:
            species = species or 'human'
            try:
                db.initialize_gene_table(species)
            except DuplicateGeneTableError:
                raise DuplicateGeneTableError(f'Gene table for {species} already exists in database. Please use a different species or set initialize_gene_table to False')
        db.conn.execute(_core_queries['insert_de_arrows'], (species,))
        db.conn.commit()
        db.conn.close()
    
    def add_experiment(
        self,
        data: object,
        metadata: MetadataInput | None = None,
        experiment_id: ExperimentId | None = None,
        **fields: ExperimentMetadataField,
    ) -> "DeArrows":
        """
        Add a new experiment to the DeArrows object. The new experiment can be provided as a DeArrow, DeArrows, or a polars or pandas DataFrame or a file path. The metadata for the new experiment can be provided as a dictionary or as keyword arguments. If no metadata is provided, the metadata from the existing DeArrows object will be used.
        The new experiment will be added to the existing DeArrows object and a new DeArrows object will be returned.
        Uses functions outside of the class.
        """
        metadata = _try_process_metadata(metadata) if metadata is not None else None
        meta = dict(self.experiment_metadata)
        if metadata is not None and fields:
            metadata.update(fields)
        if isinstance(data, DeArrow):
            df, meta, new_ids = _add_experiment_arrow(self._table.clone(), meta, self.id, data, experiment_id)
        elif isinstance(data, DeArrows):
            df, meta, new_ids = _add_experiment_arrows(self._table.clone(), meta, self.id, data, experiment_id)
        else:
            df, meta, new_ids = _add_experiment_data(self._table.clone(), meta, self.id, data, metadata, experiment_id)
        return DeArrows._from_arrow(df, meta, new_ids)
    
    def df(self) -> pl.DataFrame:
        return self._table.clone()
    
    def write_parquet(self, output_path: str, experiment_id: ExperimentId | None | list[ExperimentId] = None, compression: str = 'zstd', compression_level: int = 3) -> None:
        """
        Calls the _write_parquet function to write the DeArrows object to parquet.
        """
        experiment_id = [experiment_id] if isinstance(experiment_id, int) else experiment_id
        table = self._table.filter(pl.col('experiment_id').is_in(experiment_id)) if experiment_id is not None else self._table
        _write_parquet(table, self.experiment_metadata, output_path, compression=compression, compression_level=compression_level)

def _add_experiment_arrow(
    df: pl.DataFrame,
    meta: ExperimentMetadataMap,
    ids: list[ExperimentId],
    data: DeArrow,
    experiment_id: ExperimentId | None = None,
) -> tuple[pl.DataFrame, ExperimentMetadataMap, list[ExperimentId]]:
    frame = data._table
    if experiment_id is not None:
        frame = frame.with_columns(pl.lit(experiment_id).alias('experiment_id'))
    df.extend(frame)
    experiment_id = experiment_id or data.id
    new_ids = ids + [experiment_id]
    _check_ids(new_ids)
    new_metadata = {experiment_id: data.experiment_metadata[data.id]}
    meta.update(new_metadata)
    return (df, meta, new_ids)

def _add_experiment_arrows(
    df: pl.DataFrame,
    meta: ExperimentMetadataMap,
    ids: list[ExperimentId],
    data: DeArrows,
    experiment_id: list[ExperimentId] | None = None,
) -> tuple[pl.DataFrame, ExperimentMetadataMap, list[ExperimentId]]:
    frame = data._table
    if experiment_id is not None:
        if not isinstance(experiment_id, list):
            raise DeQuackError(f'experiment_id must be a list when adding DeArrows, not {type(experiment_id)}')
        _check_ids(ids)
        if isinstance(experiment_id, list):
            if len(experiment_id) != len(data.id):
                raise DeQuackError('Number of provided ids does not match number of existing ids')
            id_dict = dict(zip(data.id, experiment_id))
            frame = frame.with_columns(pl.col('experiment_id').replace(id_dict).alias('experiment_id'))
            new_metadata = {n: value for n, value in zip(experiment_id, data.experiment_metadata.values())}
    else:
        new_metadata = data.experiment_metadata
    df = df.extend(frame)
    new_id = experiment_id or data.id
    new_ids = ids + new_id
    _check_ids(new_ids)
    meta.update(new_metadata)
    return (df, meta, new_ids)

def _add_experiment_data(
    df: pl.DataFrame,
    meta: ExperimentMetadataMap,
    ids: list[ExperimentId],
    data: object,
    metadata: ExperimentMetadataRecord,
    experiment_id: ExperimentId,
) -> tuple[pl.DataFrame, ExperimentMetadataMap, list[ExperimentId]]:
    frame = _to_polars_table(data)
    frame = _order_columns(frame)
    if experiment_id is None:
        raise DeQuackError('id must be provided when adding a non-de_arrow table')
    if not isinstance(experiment_id, int):
        raise DeQuackError('id must be an integer when adding a non-de_arrow table')
    if metadata is None:
        raise DeQuackError('metadata must be provided when adding a non-de_arrow table')
    if not isinstance(metadata, dict):
        raise DeQuackError('metadata must be a dictionary when adding a non-de_arrow table')
    frame = frame.insert_column(0, pl.lit(experiment_id).alias('experiment_id'))
    metadata_fields, extra_info=ExperimentMetadata().to_dict(metadata, data)
    for key, value in json.loads(extra_info).items():
        metadata_fields[key]=value
    new_metadata = {experiment_id: metadata_fields}
    frame = _clean_df(frame)
    df.extend(frame)
    new_ids = ids + [experiment_id]
    meta.update(new_metadata)
    return (df, meta, new_ids)

   
    
def _write_parquet(
    df: pl.DataFrame,
    metadata: ExperimentMetadataMap,
    output_path: str,
    compression: str = 'zstd',
    compression_level: int = 3
) -> None:
    """
    Write the DeArrows object to a parquet file. The metadata will be written to a separate json file with the same name as the parquet file.
    This will be called by the write_parquet method of the DeArrows object. Will write the dataframe to a parquet file and the metadata to a json file with the same name as the parquet file.
    """
    if not isinstance(output_path, str):
        raise TypeError(f'output_path must be a string, not {type(output_path)}')
    if not output_path.endswith('.parquet'):
        raise ValueError(f'output_path must end with .parquet, not {output_path}')
    df.write_parquet(output_path, compression=compression, compression_level=compression_level)
    metadata_path = output_path.replace('.parquet', '_metadata.json')
    if len(os.path.dirname(metadata_path)) > 0:
        if not os.path.exists(os.path.dirname(metadata_path)):
            os.makedirs(os.path.dirname(metadata_path))
    with open(metadata_path, 'w') as f:
        json.dump(metadata, f, indent=4)
    logger.info(f'Wrote DeArrows object to {output_path} and metadata to {metadata_path}')



