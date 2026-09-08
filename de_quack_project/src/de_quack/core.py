"""
DeQuackling: Database manager for differential expression analysis data.

A manager for a DuckDB file used to store and query experimental and
gene result data. Key behaviors:

- Creates and manages tables (`experimental_data`, `gene_results`, `genes`).
- Provides helper methods to ingest and normalize incoming dataframes.
- Provides helper methods to query data by experiment, gene, or significance and returns wrapped DeArrow/DeArrows objects for further analysis.

Workflow summary (high level):
1. Open connection with `with DeQuackling() as db:` which creates tables/indexes.
2. Use `ingest()` to normalize and insert a DataFrame (preprocess_df).
3. Use helper methods like `get_experiment()`, `get_gene()`, `get_upregulated()`, etc. to query data.
"""
import duckdb
from duckdb import SQLExpression, CaseExpression, ColumnExpression, ConstantExpression, FunctionExpression
import re
import os
import json
import uuid
from pathlib import Path
import hashlib
import difflib
import datetime
import polars as pl
from typing import TypeAlias
from .exceptions import ProcessingError, DuplicateExperimentError, DuplicateGeneTableError, DeQuackError
from .utilities import gene_columns, ExperimentMetadata, CORE_QUERIES, _setup_logger, _try_process_metadata
import time
import importlib
core_queries = CORE_QUERIES
experiment_columns=['experiment_id', 'model', 'date', 'file', 'experiment_name', 'contrast', 'annotation_version', 'normalization', 'other_info']

ExperimentId: TypeAlias = int
ExperimentMetadataField: TypeAlias = str | int | float | bool | None | dict[str, object] | list[object]
ExperimentMetadataRecord: TypeAlias = dict[str, ExperimentMetadataField]
ExperimentMetadataMap: TypeAlias = dict[ExperimentId, ExperimentMetadataRecord]
MetadataInput: TypeAlias = dict|list[dict]|str|pl.DataFrame| None

_GENE_ALIAS_TO_COLUMN = {
    alias.lower(): canonical
    for canonical, aliases in gene_columns.items()
    for alias in aliases + [canonical]
}

logger = _setup_logger()

#make excel sheet parsing
#add more parsing compatibility for metadata

class DeQuackling:
    def __init__(self, db_path: str = 'SQL.duckdb') -> None:
        """Initialize a DuckDB-backed differential expression store.

        Args:
            db_path: DuckDB database path to open or create.
        """
        self.db_path = db_path
        self.conn = None
    
        
        
        
    def __enter__(self) -> "DeQuackling":
        """Open the DuckDB connection and create the required schema if needed.

        The core tables are `experimental_data`, `gene_results`, and `genes`.
        """
        self.conn = duckdb.connect(self.db_path)
        tables = self.conn.execute("SELECT table_name FROM information_schema.tables WHERE table_schema='main'").fetchall()
        self.tables = [table[0] for table in tables]
        if 'gene_results' not in self.tables:
            self._initialize_tables()

        self.gene_columns = self.conn.execute("SELECT * FROM information_schema.columns WHERE table_name='gene_results'").fetchall()
        self.experiment_columns = self.conn.execute("SELECT * FROM information_schema.columns WHERE table_name='experimental_data'").fetchall()
        return self
    
    
    def _initialize_tables(self) -> None:
        """Create the core tables, sequence, and indexes used by the database."""
        self.conn.execute('''CREATE SEQUENCE IF NOT EXISTS seq_experiment_id START 1''')

        # Create tables if they don't already exist.
        self.conn.execute('''CREATE TABLE IF NOT EXISTS experimental_data
                            (experiment_id INTEGER PRIMARY KEY DEFAULT nextval('seq_experiment_id'), 
                             model VARCHAR, date VARCHAR, file VARCHAR, 
                             experiment_name VARCHAR, contrast VARCHAR, annotation_version VARCHAR, normalization VARCHAR, other_info VARCHAR, data_signature VARCHAR,
                             duckdb_version VARCHAR, de_quack_version VARCHAR)''')
        
        
        
        self.conn.execute('''CREATE TABLE IF NOT EXISTS gene_results
                            (experiment_id INTEGER,
                            gene_symbol VARCHAR,
                            ensembl_id VARCHAR,
                            log2fc DOUBLE,
                            logCPM DOUBLE,
                            pvalue DOUBLE,
                            padj DOUBLE,
                            stat DOUBLE,
                            other_info VARCHAR)
                            ''')
        
        self.conn.execute('''CREATE TABLE IF NOT EXISTS genes
                            (
                            symbol VARCHAR, 
                            id VARCHAR,
                            ensembl_id VARCHAR,
                            alias_symbol VARCHAR[],
                            prev_symbol VARCHAR[], species VARCHAR,
                            PRIMARY KEY (id, species) 
                            )''')
        

        self.conn.execute("CREATE INDEX IF NOT EXISTS idx_ref_symbol ON genes (symbol)")
        self.conn.execute("CREATE INDEX IF NOT EXISTS index_ref_ensembl_id ON genes (ensembl_id)")

        
        
        
        # Create indexes
        self.conn.execute('CREATE INDEX IF NOT EXISTS idx_gene_name_lookup ON gene_results(gene_symbol, experiment_id)')
        logger.info("Initialized DuckDB schema with tables: experimental_data, gene_results, genes.")
    
    def __exit__(self, exc_type: object, exc_value: object, traceback: object) -> None:
        """Close the DuckDB connection when leaving the context manager."""
        if self.conn:
            self.conn.close()

    def fast_connect(self) -> "DeQuackling":
        """Open a connection without initializing the schema."""
        self.conn = duckdb.connect(self.db_path)
        return self
    
    def connect(self) -> "DeQuackling":
        """Ensure the database connection is open and return the manager."""
        if not self.conn:
            self.__enter__()
        return self



    def _get_data_signature(self) -> str:
        """Deterministic SHA-256 signature of the full normalized dataset,
        order-independent (rows hashed individually, then combined)."""
        view = self.conn.view('preprocessed_data')
        cols = ', '.join(f'CAST({c} AS VARCHAR)' for c in view.columns)
        result = self.conn.execute(f'''
            SELECT sha256(string_agg(row_hash, '' ORDER BY row_hash))
            FROM (
                SELECT sha256(concat_ws('\x1f', {cols})) AS row_hash
                FROM view
            )
        ''').fetchone()
        logger.info('Computed data signature for normalized dataset.')
        return result[0]

   
    
        

    def _create_temp_view(self, columns: dict[str, str] | None = None) -> duckdb.DuckDBPyRelation:
        """Create a normalized temporary view from the registered input data.
            Uses a mix of provided column mappings (if any) and the default gene columns found in utilities.gene_columns. Raises ProcessingError if no recognizable gene columns are found.
            Returns the view as a relation to use in the insert query. On the event that basemean is detected, it will be used to compute logCPM if logCPM is not present.
            Basemean will then be added to the view as a column, but will be in the extra columns section of the other_info JSON column.
            If an ambiguous gene name is detected (eg. a column that could be symbol or Ensembl ID), it will be treated as a symbol unless the data is detected to be Ensembl IDs (eg. starts with ENS).
            """
        if columns is None:
            columns = {}
        temp_view = self.conn.table('preprocessed_data')
        columns_info = [(column_name, _GENE_ALIAS_TO_COLUMN.get(column_name.lower())) for column_name in temp_view.columns]
        for key, value in columns.items():
            if key in temp_view.columns and value in gene_columns:
                if (key, value) not in columns_info:
                    columns_info.append((key, value))
            else:
                raise ProcessingError(f"Invalid column mapping: {key} -> {value}. Ensure the column exists in the data and the value is a valid gene column alias.")


        if not columns_info:
            raise ProcessingError("The registered table contains no columns.")

        has_gene_symbol = any(real_col == 'gene_symbol' for _, real_col in columns_info)

        mapped_incoming = []
        all_columns = []
        mapped_columns = []
        expressions=[]

        for incoming_col, real_col in columns_info:
            all_columns.append(incoming_col)
            if real_col == 'gene_name':
                sample = temp_view.select(ColumnExpression(incoming_col)).limit(1).fetchone()
                sample_val = sample[0] if sample is not None else None
                if sample_val is not None and re.match(r'^ENS', str(sample_val)):
                    mapped_incoming.append(incoming_col)
                    expressions.append(ColumnExpression(incoming_col).alias('ensembl_id'))
                    mapped_columns.append('ensembl_id')
                elif not has_gene_symbol:
                    mapped_incoming.append(incoming_col)
                    expressions.append(ColumnExpression(incoming_col).alias('gene_symbol'))
                    mapped_columns.append('gene_symbol')
                else:
                    pass
            elif real_col is not None:
                if real_col != 'base_mean':
                    mapped_incoming.append(incoming_col)
                expressions.append(ColumnExpression(incoming_col).alias(real_col))
                mapped_columns.append(real_col)

        if not mapped_incoming:
            raise ProcessingError(
                'No recognizable gene result columns were found in the input. '
                'Verify the header names against gene_columns aliases.'
            )
        for col in ['gene_symbol', 'ensembl_id', 'log2fc', 'padj', 'pvalue', 'stat', 'logCPM']:
            if col not in mapped_columns:
                expressions.append(ConstantExpression(None).alias(col))
                mapped_columns.append(col)
    
        extra_columns = [col for col in all_columns if col not in mapped_incoming]
        if extra_columns:
            struct_items=[]
            for column in extra_columns:
                struct_items.append(ConstantExpression(column))
                struct_items.append(ColumnExpression(column))
            json_expr = FunctionExpression('json_object', *struct_items)
            expressions.append(json_expr.alias('other_info'))
        else:
            expressions.append(ConstantExpression(None).alias('other_info'))
        
        
        
        if 'base_mean' not in mapped_columns:
            expressions.append(SQLExpression('NULL').alias('base_mean'))

        if not mapped_incoming:
            raise ProcessingError(
                'No recognizable gene result columns were found in the input. '
                'Verify the header names against gene_columns aliases.'
            )
        temp_view=temp_view.select(*expressions)
        return temp_view




    def _preprocess(self, info: object) -> None:
        """Register an input file or table-like object as `preprocessed_data`.
            Supports CSV, TSV, Parquet, DuckDB relation, Polars DataFrame, Pandas DataFrame, DeArrow, or DeArrows.
        """

        if isinstance(info, str) and os.path.isfile(os.path.abspath(info)):
            info=os.path.abspath(info)
            if info.lower().endswith('.parquet'):
                self.conn.execute('DROP VIEW IF EXISTS preprocessed_data')
                self.conn.read_parquet(info).create_view('preprocessed_data')
                return
            elif info.lower().endswith('xls') or info.lower().endswith('xlsx'):
                try:
                    import openpyxl
                except ImportError:
                    raise ImportError('openpyxl is required to read Excel files. Please install it with `pip install openpyxl`.')
                self.conn.execute('DROP VIEW IF EXISTS preprocessed_data')
                df = pl.read_excel(info, sheet_name=0)
                self.conn.register('info', df.to_arrow())
                self.conn.execute('CREATE TEMP VIEW preprocessed_data AS SELECT * FROM info')
                return
            else:
                self.conn.execute('DROP VIEW IF EXISTS preprocessed_data')
                self.conn.from_csv_auto(info, header=True).create_view("preprocessed_data")
                return
        if isinstance(info, str) and not os.path.isfile(os.path.abspath(info)):
            raise ProcessingError(f"{os.path.abspath(info)} is not a valid file path.")
        

        self.conn.execute('DROP VIEW IF EXISTS preprocessed_data')
        from .arrow import DeArrow, DeArrows
        if isinstance(info, DeArrow) or isinstance(info, DeArrows):
            self.conn.register('info', info._table)
            self.conn.execute('CREATE TEMP VIEW preprocessed_data AS SELECT * FROM info')
            return
        
        if isinstance(info, duckdb.DuckDBPyRelation):
            self.conn.register('info', info)
            self.conn.execute('CREATE TEMP VIEW preprocessed_data AS SELECT * FROM info')
            return

        try:
            import pandas as pd
            if isinstance(info, pd.DataFrame):
                self.conn.register('info', info)
                self.conn.execute('CREATE TEMP VIEW preprocessed_data AS SELECT * FROM info')
                return
        except ImportError:
            pass

        try:
            import pyarrow as pa
            if isinstance(info, pa.Table):
                self.conn.register('info', info)
                self.conn.execute('CREATE TEMP VIEW preprocessed_data AS SELECT * FROM info')
                return
        except ImportError:
            pass

        raise ProcessingError(f'type {type(info)} is not supported. Provide a pandas DataFrame, DuckDB relation, de_arrow object, or path to a CSV/TSV/Parquet file.')
    
    def initialize_gene_table(self, species: str = 'human') -> None:
        """Load the reference gene table for `human` or `mouse`.

        `human` loads HGNC data; `mouse` loads MGI data.
        """
        if species.lower() == 'human':
            self._insert_human_genes()
        elif species.lower() == 'mouse':
            self._insert_mouse_genes()
        else:
            raise ProcessingError(f"Species '{species}' is not supported. Only 'human' and 'mouse' are currently supported.")
    
    def _insert_mouse_genes(self) -> None:
        """Load mouse reference gene annotations into the `genes` table.

        The data comes from the bundled mouse gene mapping query.
        """
        from .arrow import _get_gene_table
        table_path = _get_gene_table('mouse')
        
       
        try:
            self.conn.execute('INSERT INTO genes SELECT * FROM read_parquet(?)', (str(table_path),))
            self.conn.commit()
        except duckdb.ConstraintException as e:
            raise DuplicateGeneTableError(
                'Gene reference data appears to have already been loaded or the gene table has duplicate entries.'
            ) from e
        logger.info('Mouse gene reference data loaded into the `genes` table.')

    

    def _insert_human_genes(self) -> None:
        """Load human reference gene annotations into the `genes` table.

        The data comes from the bundled human gene mapping query.
        """
        from .arrow import _get_gene_table
        table_path = _get_gene_table('human')
        try:
            self.conn.execute('INSERT INTO genes SELECT * FROM read_parquet(?)', (str(table_path),))
            self.conn.commit()

        except duckdb.ConstraintException as e:
            raise DuplicateGeneTableError(
                'Gene reference data appears to have already been loaded or the gene table has duplicate entries.'
            ) from e
        logger.info('Human gene reference data loaded into the `genes` table.')

    def _create_experiment(
        self,
        metadata_fields: ExperimentMetadataRecord,
        data_signature: str,
        other_info: str,
    ) -> ExperimentId:
        """Insert one experiment row and return the generated experiment id."""
        duplicate_ids = self.conn.execute(
            "SELECT experiment_id FROM experimental_data WHERE data_signature = ?", 
            (data_signature,)
        ).fetchall()
        if len(duplicate_ids) > 0:
            raise DuplicateExperimentError(
                f'Data is identical to data in experiment id: {duplicate_ids[0][0]}. Duplicate experiments are not allowed.'
            )
        result = self.conn.execute(
            "INSERT INTO experimental_data (model, date, file, experiment_name, contrast, annotation_version, normalization, other_info, data_signature, duckdb_version, de_quack_version) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?) RETURNING experiment_id",
            (metadata_fields.get('model'), metadata_fields.get('date'), metadata_fields.get('file'), metadata_fields.get('experiment_name'), metadata_fields.get('contrast'), metadata_fields.get('annotation_version'), metadata_fields.get('normalization'), other_info, data_signature, importlib.metadata.version('duckdb'), importlib.metadata.version('de_quack'))
        ).fetchall()
        self.conn.commit()
        
        return result[0][0]
    
    def write_parquet(self, output_path: str, experiment_id: ExperimentId | None | list[ExperimentId] = None) -> None:
        """Write gene results and experiment metadata to a Parquet file for external use.
        Args:
            output_path: Path to write the Parquet file.
            experiment_id: Optional experiment_id to filter the results. If None, all experiments are written.
        """
        logger.info(f"Writing gene results and experiment metadata to {output_path}.")

        if len(os.path.dirname(output_path)) > 0:
            os.makedirs(os.path.dirname(output_path), exist_ok=True)
        output_path = os.path.abspath(output_path)
        if experiment_id is None:
            table = self.conn.table('gene_results')
            table.write_parquet(output_path)
            metadata = self.conn.execute(f"SELECT * FROM experimental_data").fetchall()
            experiment_id = [n for n in range(1, len(metadata)+1)]
        elif isinstance(experiment_id, list):
            try:
                for i in range(len(experiment_id)):
                    experiment_id[i] = int(experiment_id[i])
            except ValueError:
                raise ValueError("All experiment IDs must be integers.")
            table = self.conn.table('gene_results').filter(ColumnExpression('experiment_id').isin(*experiment_id))
            table.write_parquet(output_path)
            metadata = self.conn.execute(f"SELECT * FROM experimental_data WHERE experiment_id = ANY($1)", (experiment_id,)).fetchall()
        else:
            try:
                experiment_id = int(experiment_id)
            except ValueError:
                raise ValueError("Experiment ID must be an integer.")
            table = self.conn.table('gene_results').filter(ColumnExpression('experiment_id') == experiment_id)
            table.write_parquet(output_path)
            metadata = self.conn.execute(f"SELECT * FROM experimental_data WHERE experiment_id = ?", (experiment_id,)).fetchall()
        meta_list = [dict(zip(experiment_columns, row)) for row in metadata]
        meta_dict = _to_metadata(meta_list)
        meta_json = json.dumps(meta_dict, indent=4)
        meta_path = os.path.splitext(output_path)[0] + '_metadata.json'
        with open(meta_path, 'w') as f:
            f.write(meta_json)
        logger.info(f'Gene results written to {output_path}.')
        logger.info(f'Experiment metadata written to {meta_path}.')

    def _insert_gene_results(self, experiment_id: ExperimentId, view: duckdb.DuckDBPyRelation, species: str | None = None) -> None:
        """Insert normalized gene-result rows for one experiment. Uses second query to get unmatched genes and log a warning if any are found."""



        self.conn.execute('''
        INSERT INTO gene_results (experiment_id, gene_symbol, ensembl_id, log2fc, logCPM, pvalue, padj, stat, other_info)
        SELECT $1 AS experiment_id,
            COALESCE(ens.symbol, sym.symbol, e.gene_symbol) AS gene_symbol,
            COALESCE(ens.ensembl_id, sym.ensembl_id, e.ensembl_id) AS ensembl_id,
            COALESCE(TRY_CAST(e.log2fc AS DOUBLE), NULL) AS log2fc,
            COALESCE(TRY_CAST(e.logCPM AS DOUBLE), log2(e.base_mean + 1), NULL) AS logCPM,
            COALESCE(TRY_CAST(e.pvalue AS DOUBLE), NULL) AS pvalue,
            COALESCE(TRY_CAST(e.padj AS DOUBLE), NULL) AS padj,
            COALESCE(TRY_CAST(e.stat AS DOUBLE), NULL) AS stat,
            e.other_info
        FROM view as e
        LEFT JOIN genes ens ON  (
            e.ensembl_id IS NOT NULL AND ens.ensembl_id = UPPER(TRIM(e.ensembl_id))
        )

        LEFT JOIN genes sym ON 
            (e.gene_symbol IS NOT NULL AND sym.symbol = e.gene_symbol) OR
            e.gene_symbol IS NOT NULL AND sym.prev_symbol IS NOT NULL AND list_contains(sym.prev_symbol, e.gene_symbol)
        
        
        WHERE (sym.species IS NULL OR sym.species = $2) AND (ens.species IS NULL OR ens.species = $2)
        ''', (experiment_id, species))

        mismatch_count = self.conn.execute("""
        SELECT COUNT(*)
        FROM view AS e
        LEFT JOIN genes ens ON (
            e.ensembl_id IS NOT NULL AND ens.ensembl_id = UPPER(TRIM(e.ensembl_id))
        )
        LEFT JOIN genes sym ON (
            e.gene_symbol IS NOT NULL AND sym.symbol = e.gene_symbol
        ) OR (
        e.gene_symbol IS NOT NULL AND sym.prev_symbol IS NOT NULL AND list_contains(sym.prev_symbol, e.gene_symbol)
        )
        WHERE (sym.species IS NULL OR sym.species = $1) AND (ens.species IS NULL OR ens.species = $1)
        AND ens.ensembl_id IS NULL AND sym.symbol IS NULL
        """, (species,)).fetchone()[0]

        if mismatch_count > 0:
            logger.warning(f"{mismatch_count} gene results could not be matched to the reference gene table for species '{species}'. Check the 'gene_symbol' and 'ensembl_id' columns for potential mismatches.")

        self.conn.commit()
        

    def ingest(
        self,
        info: object,
        metadata: MetadataInput | None = None,
        species: str = 'human',
        columns: dict[str, str] | None = None,
        **kwargs: ExperimentMetadataField,
    ) -> None:
        """Normalize and ingest differential expression results into DuckDB.

        Args:
            metadata: Experiment metadata for the ingested table. Can be a dictionary, a string path to a JSON or CSV file, or a polars or pandas DataFrame.
            columns: Optional explicit column remapping for gene fields.
            **kwargs: Additional metadata fields merged into `metadata`.
        """
        _try_process_metadata(metadata)
        if isinstance(metadata, list) and len(metadata) > 1:
            logger.warning('Multiple metadata records provided. Only the first will be used for this insertion.')
        metadata = metadata[0] if isinstance(metadata, list) else metadata
        if columns is None:
            columns = {}
        if metadata is None:
            metadata = {}
        if kwargs:
            # Check for typos in the kwargs keys
            parameters_and_keys = experiment_columns + ['metadata, columns, species']
            for key in kwargs.items():
                if difflib.get_close_matches(key, parameters_and_keys, n=1, cutoff=0.8):
                    raise ProcessingError(f"Invalid metadata field: {key}. Did you mean '{difflib.get_close_matches(key, parameters_and_keys, n=1)[0]}'?")
            metadata.update(kwargs)
        metadata_fields, other_info = ExperimentMetadata().to_dict(metadata, info)
        self._preprocess(info) 

        data_signature = self._get_data_signature()

        experiment_id = self._create_experiment(metadata_fields, data_signature, other_info)

        view=self._create_temp_view(columns)

        self._insert_gene_results(experiment_id, view, species)
        
    def get_experiment(
        self,
        experiment_id: ExperimentId | None = None,
        name: str | None = None,
        model: str | None = None,
        annotation_version: str | None = None,
        normalization: str | None = None,
        date: str | None = None,
        contrast: str | None = None,
        file: str | None = None,
    ) -> "DeArrow | DeArrows":
        """Return experiments matching the provided metadata filters."""
        try:
            date = datetime.datetime.strptime(date, '%Y-%m-%d').date() if date else None
        except ValueError:
            raise ProcessingError(f"Invalid date format: {date}. Expected format is YYYY-MM-DD.")
        rel = self.conn.execute(core_queries['get_experiment'], [experiment_id, date, model, file, name, contrast, annotation_version, normalization]).fetchall()
        metadata_list = [dict(zip(experiment_columns, row)) for row in rel]
        metadata_fields, ids = _to_metadata(metadata_list)
        self.conn.execute(f'CREATE OR REPLACE TEMP TABLE ids AS SELECT UNNEST({ids}) AS id')
        df = self.conn.execute('SELECT * exclude(ids.id) FROM gene_results g JOIN ids ON g.experiment_id = ids.id').pl()
        return _get_de_arrows(df, metadata_fields, ids)
        
    
    def get_significant_genes(self, log2fc: float = 1, padj: float = 0.05, logCPM: float = 1, experiment_id: ExperimentId | None = None) -> "DeArrow | DeArrows":
        """Return genes meeting the log2 fold change, padj, and expression cutoffs."""
        df = self.conn.execute('SELECT * FROM gene_results WHERE abs(log2fc) > $1 AND padj < $2 AND logCPM > $3 AND (experiment_id = $4 OR $4 IS NULL)', [log2fc, padj, logCPM, experiment_id]).pl()
        return self._polars_to_de_arrows(df)

    def get_upregulated(self, log2fc: float = 1, padj: float = 0.05, logCPM: float = 1, experiment_id: ExperimentId | None = None) -> "DeArrow | DeArrows":
        """Return genes with positive log2 fold change above the threshold."""
        df = self.conn.execute('SELECT * FROM gene_results WHERE log2fc > $1 AND padj < $2 AND logCPM > $3 AND (experiment_id = $4 OR $4 IS NULL)', [log2fc, padj, logCPM, experiment_id]).pl()
        arrow = self._polars_to_de_arrows(df)
        return arrow
        
    
    def get_downregulated(self, log2fc: float = -1, padj: float = 0.05, logCPM: float = 1, experiment_id: ExperimentId | None = None) -> "DeArrow | DeArrows":
        """Return genes with negative log2 fold change below the threshold."""
        df = self.conn.execute('SELECT * FROM gene_results WHERE log2fc < $1 AND padj < $2 AND logCPM > $3 AND (experiment_id = $4 OR $4 IS NULL)', [log2fc, padj, logCPM, experiment_id]).pl()
        return self._polars_to_de_arrows(df)

    def get_gene(self, gene_symbol: str | None = None, ensembl_id: str | None = None, experiment_id: ExperimentId | None = None) -> "DeArrow | DeArrows":
        """Return genes matching the provided symbol, Ensembl id, and/or experiment id."""
        rel = self.conn.execute('SELECT * FROM gene_results WHERE (gene_symbol = $1 OR $1 IS NULL) AND (ensembl_id = $2 OR $2 IS NULL) AND (experiment_id = $3 OR $3 IS NULL)', [gene_symbol, ensembl_id, experiment_id]).pl()
        return self._polars_to_de_arrows(rel)

    def delete_experiment(
        self,
        experiment_id: ExperimentId | None = None,
        name: str | None = None,
        model: str | None = None,
        annotation_version: str | None = None,
        normalization: str | None = None,
        date: str | None = None,
        contrast: str | None = None,
        file: str | None = None,
    ) -> None:
        """Delete experiments and their gene-result rows for the matching filters."""
        try:
            date = datetime.datetime.strptime(date, '%Y-%m-%d').date() if date else None
        except ValueError:
            raise ProcessingError(f"Invalid date format: {date}. Expected format is YYYY-MM-DD.")
        try:
            self.conn.begin()
            ids = self.conn.execute(core_queries['find_delete_experiment'], (experiment_id, date, model, file, name, contrast, annotation_version, normalization)).fetchall()
            ids = [row[0] for row in ids]
            deleted_genes = self.conn.execute('SELECT COUNT(*) FROM gene_results WHERE experiment_id = ANY(?)', (ids,)).fetchone()[0]
            self.conn.execute(core_queries['delete_gene_results'], (ids,))
            self.conn.execute(core_queries['delete_experiment'], (ids,))
            logger.info(f'Deleted {deleted_genes} gene result rows and {len(ids)} experiment metadata rows for experiment_id(s): {ids}.')
            logger.info(f'Deleted {len(ids)} experiment(s) with experiment_id(s): {ids}.')
            self.conn.commit()
        except Exception as e:
            self.conn.rollback()
            raise DeQuackError(f"Failed to delete experiment(s): {e}") from e
            
    
    def _polars_to_de_arrows(self, df: pl.DataFrame) -> "DeArrow | DeArrows":
        """Convert a Polars result frame into a `DeArrow` or `DeArrows` object."""
        if len(df.columns) == 0:
            return self._get_de_arrows(df, {}, [])
        ids = df['experiment_id'].unique().to_list()
        meta_rel = self.conn.execute(
            '''
            SELECT
                experiment_id,
                model,
                date,
                file,
                experiment_name,
                contrast,
                annotation_version,
                normalization,
                other_info
            FROM experimental_data
            WHERE experiment_id = ANY(?)
            ''',
            [ids]
        )
        metadata_list = [dict(zip(experiment_columns, row)) for row in meta_rel.fetchall()]
        metadata_fields, ids = _to_metadata(metadata_list)
        return _get_de_arrows(df, metadata_fields, ids)
        
            



    
        
        

    

    def close(self) -> None:
        """Close the active DuckDB connection and clear the handle."""
        if self.conn:
            self.conn.close()
            self.conn = None



    
def _get_de_arrows(arrows: object, metadata: ExperimentMetadataMap, ids: list[ExperimentId]) -> "DeArrow | DeArrows":
    """ Imports arrow at the function level to avoid circular import issues. Returns a DeArrow or DeArrows object based on the number of experiment ids."""
    from .arrow import DeArrow, DeArrows
    if len(ids) < 2:
        return DeArrow._from_arrow(arrows, metadata, ids)
    return DeArrows._from_arrow(arrows, metadata, ids)

def _to_metadata(metadata_list: list[ExperimentMetadataRecord]) -> tuple[ExperimentMetadataMap, list[ExperimentId]]:
    """Flatten experiment metadata rows and merge any JSON `extra_info` payloads."""
    for meta in metadata_list:
        if meta.get('extra_info'):
            try:
                extra = json.loads(meta['extra_info'])
                for key, value in extra.items():
                    meta[key] = value
            except json.JSONDecodeError:
                raise ProcessingError (f"Failed to parse 'extra_info' JSON for experiment_id {meta.get('experiment_id')}. Ensure it is valid JSON.")
        meta.pop('extra_info', None)
    ids = [meta.pop('experiment_id') for meta in metadata_list]
    metadata_fields = {exp_id: meta for meta, exp_id in zip(metadata_list, ids)}
    return metadata_fields, ids





    

        