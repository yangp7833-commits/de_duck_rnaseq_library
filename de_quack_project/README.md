# de_quack

de_quack is a Python library for storing, normalizing, and querying differential expression (DE) results with DuckDB. It supports:

- schema-backed storage for experiment metadata and gene-level DE results
- automatic gene symbol / Ensembl ID normalization
- species-aware gene annotation via bundled human and mouse reference tables
- in-memory and file-based data workflows with Polars/DataFrame support
- volcano plot generation for DE result exploration

Purpose:

- while there are many database managers for omics data, there are few that are specifically designed for differential expression results or independent work
- current tools are difficult to start up and configurate, and often require a large database server or cloud service, leading many to use custom file workflows that are not easily shared or reproducible
- de_quack aims to fill this gap by providing a simple and integrative tool for DE result storage and analysis, with a focus on reproducibility and ease of use, while keeping the consistency of traditional database managers through column name checking and schema enforcement


## Installation

```bash
pip install de_quack
```

Optional plotting extras:

```bash
pip install "de_quack[plotting]"
```
## Changes from V1

- The library has been refactored to use Polars as the main data engine, replacing Pandas. This change improves performance and memory efficiency, especially for large datasets.
- The previous reticulate bridge has been removed, and the library is now fully implemented in Python. This simplifies the codebase and reduces dependencies.
- The `DeArrow` and `DeArrows` classes have been introduced to represent single and multiple experiments, respectively. These classes provide a more intuitive interface for working with DE results. They are also built on the Polars Dataframe, which allows for efficient data manipulation and querying.
- The `DeQuackling` class has been updated to work with the new Polars-based data model. It now provides methods for ingesting, querying, and managing DE results in a DuckDB database.
- The `volcano_plot` function has been updated to work with Polars DataFrames, and it now supports additional customization options for plotting.
- The library now includes improved error handling and input validation, ensuring that users receive informative messages when issues arise.

## Supported inputs

The library accepts many table-like objects and file paths:

- CSV / TSV / Parquet files
- pandas DataFrames
- Polars DataFrames
- Polars LazyFrames
- DuckDB relations
- Excel files (xlsx) 
- `DeArrow` and `DeArrows` objects

With version 2.1, the library now supports CSV, TSV, and JSON files for metadata input along with Polars and pandas Dataframes. 

Note: Although these input types are supported, they are detected by first trying an import of the relevant library and then checking the type of the object. Therefore, even though Pandas is supported, it is not a dependency.
## Core data model

de_quack stores results in a DuckDB database with three main tables:

- `experimental_data`: metadata for each DE experiment
- `gene_results`: gene-level statistics and annotations
- `genes`: reference gene symbol / Ensembl mapping for a species

The canonical gene columns used by the library are:

- `gene_symbol`
- `ensembl_id`
- `log2fc`
- `logCPM`
- `pvalue`
- `padj`
- `stat`
- `other_info`

Common aliases are recognized automatically, such as `symbol`, `gene_id`, `logfc`, `p_adj`, `qvalue`, and `basemean`. Here are the official recognized aliases for each column, note that base_mean and gene_name are not canonical columns and are changed to fit the canonical schema (read below for more details):
    `gene_symbol`: gene_symbol, symbol, hgnc_symbol, genesymbol
    `gene_name`: gene, gene_name, genename, name
    `ensembl_id`: ensembl_id, ensembl, gene_id, geneid, ensembl_gene_id, target_id, feature_id
    `pvalue`: pvalue, p_value, p.value, p-value, pval, p.val, p_val
    `padj`: padj, p.adj, p_adj, p.adjusted, p_adjusted, qvalue, q_value, q.val, q_val, fdr
    `stat`: stat, f, lr, t
    `log2fc`: log2foldchange, log2fc, log2_fc, log2.fc, logfc, log_fc, log.fc
    `base_mean`: basemean, base_mean,  aveexpr, tpm, fpkm
    `logCPM`: logcpm, log_cpm, log.cpm

If base_mean is detected, the logCPM is calculated from it using logCPM = log2(base_mean + 1), and the original base_mean data is stored in the other_info column. The gene_name column is used to detect ambiguous column names. If this column is detected, a regex will be performed to determine if it is a gene symbol or an ensembl id. If it is a gene symbol, it will be stored in the gene_symbol column, and if it is an ensembl id, it will be stored in the ensembl_id column. The original gene_name data is stored in the other_info column.

## Quick start

### Python

```python
from de_quack import DeQuackling, volcano_plot

with DeQuackling("results.duckdb") as db:
    db.initialize_gene_table("human")

    db.ingest(
        "data.txt",
        metadata={
            "experiment_name": "Test experiment",
            "contrast": "treatment_vs_control",
            "model": "DESeq2",
            "normalization": "voom",
        },
        species="human",
        columns={"log2FoldChange": "log2fc"},
    )

    genes = db.get_significant_genes(log2fc=1, padj=0.05, logCPM=1)
    print(genes.df().head())

    up = db.get_upregulated(log2fc=1, padj=0.05, logCPM=1)
    print(up.df().head())

    one_gene = db.get_gene(gene_symbol="AKT1")
    print(one_gene.df())

    volcano_plot(genes, padj=0.05, log2fc=1, title="Example volcano plot", show=True)
```

## Important classes and functions

### `DeQuackling`

The main database manager. It opens a DuckDB file, creates the necessary schema, stores DE experiments, and offers query helpers.

`DeQuackling(db_path: str = "SQL.duckdb")`

Methods:

- `__enter__()` / `__exit__()`: open and close the DuckDB connection via context manager usage.
- `connect()`: ensure the database is connected.
- `close()`: close the active DuckDB connection.
- `initialize_gene_table(species: str = "human")`: load the bundled gene reference table for `human` or `mouse`.
- `ingest(info, metadata=None, species="human", columns=None, **kwargs)`: normalize and insert a DE result table into the database.
- `get_experiment(experiment_id=None, name=None, model=None, annotation_version=None, normalization=None, date=None, contrast=None, file=None)`: fetch matching experiment metadata and rows.
- `get_significant_genes(log2fc=1, padj=0.05, logCPM=1, experiment_id=None)`: return genes passing combined fold-change, adjusted-p-value, and expression cutoffs.
- `get_upregulated(log2fc=1, padj=0.05, logCPM=1, experiment_id=None)`: return genes with positive log2 fold change above threshold.
- `get_downregulated(log2fc=-1, padj=0.05, logCPM=1, experiment_id=None)`: return genes with negative log2 fold change below threshold.
- `get_gene(gene_symbol=None, ensembl_id=None, experiment_id=None)`: fetch genes by symbol, Ensembl ID, or experiment.
- `delete_experiment(experiment_id=None, name=None, model=None, annotation_version=None, normalization=None, date=None, contrast=None, file=None)`: delete matching experiments and their gene rows.
- `write_parquet(output_path, experiment_id=None)`: write gene results and experiment metadata to a Parquet file plus `_metadata.json`.

Example:

```python
with DeQuackling("results.duckdb") as db:
    db.initialize_gene_table("human")
    db.ingest(
        "data.txt",
        metadata={"experiment_name": "RNA-seq experiment", "contrast": "A_vs_B"},
        species="human",
        columns={"log2FoldChange": "log2fc"},
    )

    dat = db.get_experiment(name="RNA-seq experiment")
    print(dat.df())
```

### `DeArrow`

Represents a single experiment as a Polars-backed table with an attached experiment metadata map.

`DeArrow(info, experiment_id=None, metadata=None, heal_genes=False, species=None, columns=None, ignore_errors=False, **fields)`

Important methods:

- `get_significant_genes(log2fc=1, padj=0.05, logCPM=1)`: filter significant genes for the experiment.
- `get_downregulated(log2fc=-1, padj=0.05, logCPM=1)`: filter downregulated genes.
- `get_upregulated(log2fc=1, padj=0.05, logCPM=1)`: filter upregulated genes.
- `get_gene(gene_symbol=None, ensembl_id=None, experiment_id=None)`: fetch gene rows by ID.
- `set_id(experiment_id)`: assign or change the experiment ID while updating metadata.
- `add_experiment(data, metadata=None, experiment_id=None, **fields)`: add another experiment and return a `DeArrows` object.
- `insert(file, initialize_gene_table=False, species=None)`: insert the experiment into a DuckDB database file.
- `df()`: return a cloned Polars DataFrame.
- `write_parquet(output_path, compression="zstd", compression_level=3)`: write the experiment to Parquet.

Example:

```python
from de_quack import DeArrow

arrow = DeArrow(
    "data.txt",
    metadata={"experiment_name": "Treatment", "contrast": "drug_vs_control"},
    heal_genes=True,
    species="human",
)

sig = arrow.get_significant_genes(log2fc=1, padj=0.05, logCPM=1)
print(sig.df().head())
```

### `DeArrows`

Represents multiple experiments as one combined Polars-backed table while maintaining per-experiment metadata.

`DeArrows(*args, columns=None, metadata=None, ids=None, keep_ids=False, heal_genes=False, species="human", ignore_errors=False)`

Important methods:

- `get_experiment(experiment_id=None, name=None, model=None, annotation_version=None, normalization=None, date=None, contrast=None, file=None)`: filter by experiment metadata.
- `get_significant_genes(log2fc=1, padj=0.05, logCPM=1)`: filter across all experiments.
- `get_downregulated(log2fc=-1, padj=0.05, logCPM=1)`: filter downregulated genes.
- `get_upregulated(log2fc=1, padj=0.05, logCPM=1)`: filter upregulated genes.
- `get_gene(gene_symbol=None, ensembl_id=None, experiment_id=None)`: select rows by gene symbol or Ensembl ID.
- `set_id(ids)`: replace experiment IDs using either a list or a mapping.
- `add_experiment(data, metadata=None, experiment_id=None, **fields)`: append a new experiment to the collection.
- `df()`: access a clone of the underlying DataFrame.
- `write_parquet(output_path, experiment_id=None, compression="zstd", compression_level=3)`: write one or more experiments to Parquet.

Example:

```python
from de_quack import DeArrows

arrows = DeArrows(
    "data.txt",
    "data1.txt",
    metadata=[
        {"experiment_name": "Experiment A"},
        {"experiment_name": "Experiment B"},
    ],
    heal_genes=True,
    species="human",
)

print(arrows.get_upregulated(log2fc=1, padj=0.05, logCPM=1).df().head())
```

### `volcano_plot`

Create a volcano plot from a `DeArrow`, `DeArrows`, Polars DataFrame, or pandas DataFrame.

`volcano_plot(df, padj=0.05, log2fc=1, title=None, show=False, label_genes=0, label_type="ensembl_id", insignificant_color="grey", upregulated_color="red", downregulated_color="blue", file=None, **labeling_kwargs)`

Plot behavior:

- uses `log2fc` and `padj` columns
- highlights upregulated genes in red, downregulated genes in blue, and insignificant genes in grey
- supports optional text labels for the strongest genes

Example:

```python
from de_quack import volcano_plot

fig = volcano_plot(
    genes,
    padj=0.05,
    log2fc=1,
    title="Differential expression",
    label_genes=10,
    label_type="gene_symbol",
    file="volcano.png",
)
```

## Metadata handling

Metadata is passed as a dictionary to `ingest`, `DeArrow`, and `DeArrows`. The library normalizes common field names, including:

- `experiment_name`
- `date`
- `contrast`
- `file`
- `normalization`
- `model`
- `annotation_version`

Any additional metadata keys are stored in the `other_info` / JSON payload alongside the core fields.

## Query conventions

The library uses explicit helper methods instead of a generic SQL-like query DSL. The most useful accessors are:

- `db.get_experiment(...)`
- `db.get_gene(...)`
- `db.get_significant_genes(...)`
- `db.get_upregulated(...)`
- `db.get_downregulated(...)`

For `DeArrow` and `DeArrows`, the same filtering patterns are available on the object itself.

## Example workflow

```python
from de_quack import DeQuackling, volcano_plot

with DeQuackling("results.duckdb") as db:
    db.initialize_gene_table("human")
    db.ingest(
        "data.txt",
        metadata={
            "experiment_name": "condition_a_vs_b",
            "contrast": "A_vs_B",
            "model": "DESeq2",
            "normalization": "TMM",
        },
        species="human",
    )

    df = db.get_significant_genes(log2fc=1, padj=0.05, logCPM=1).df()
    print(df.head())

    volcano_plot(df, padj=0.05, log2fc=1, label_genes=10, label_type="gene_symbol", file="volcano.png")
```

## Notes

- `DeQuackling` enforces duplicate experiment protection using a deterministic data signature.
- `initialize_gene_table("human")` and `initialize_gene_table("mouse")` load bundled reference tables included with the package.
- The code is built around DuckDB and Polars, so it is efficient for larger DE datasets.

## Scope/Limitations
- The library is designed for differential expression analysis, not other omics data (although support for other data types could be added in the future).
- The gene table mapping is currently limited to human and mouse, but more species will be added in future releases.
- The library uses DuckDB, an embedded SQL engine. Therefore it is suited only to local analysis and independent research, not large scale multi-user management.

## Testing

- This library uses `pytest` for testing. The tests and the relevant test datasets are included in the tests directory.
- These tests are split into 3 test files:
    - `core_test.py`: tests for the DeQuackling class, covering ingestion of csv and parquet files, querying, deleting experiments, and writing to parquet
    - `arrow_test.py`: tests for the DeArrow/DeArrows classes, covering ingestion of csv and parquet files, querying, adding experiments, and writing to parquet
    - `plot_test.py`: tests for volcano_plot function, specifically plotting, labeling, and saving the plot to a file
- To run the tests, use the following command:

```bash
pytest -v tests/ #specify the file if you wish to run only one test
```
## License

MIT

