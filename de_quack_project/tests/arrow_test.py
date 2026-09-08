from de_quack import DeArrow, DeArrows, DeQuackling
import pytest
import os

class TestArrow():
    def test_arrow_creation(self):
        arrow_human = DeArrow('data.txt', metadata = {'experiment_name': 'Test Experiment'}, heal_genes = True, species = 'human')
        arrow_mouse = DeArrow('data1.txt', metadata = {'experiment_name': 'Test Experiment 2'})
        arrow_human_id = arrow_human.set_id(4)
        assert arrow_human._table.height > 0 and arrow_mouse._table.height > 0
        assert arrow_human.get_gene('AKT1').height > 0 and arrow_mouse.get_gene(ensembl_id = 'ENSMUSG00000000017').height > 0
        assert arrow_human_id.id == 4
    
    
    def test_arrows_creation(self):
        arrows = DeArrows('data.txt', 'data1.txt', heal_genes = True, species = 'human', metadata = [{'experiment_name1': 'Test Experiment'}, {'experiment_name': 'Test Experiment 2'}])
        arrows_with_ids = DeArrows(arrows, 'data.txt', metadata = {'experiment_name': 'Test Experiment 2'}, ids = [3], keep_ids = True)
        assert arrows._table.height > 0
        assert arrows.get_gene('AKT1').height > 0 and arrows.get_gene(ensembl_id = 'ENSMUSG00000000017').height > 0
        assert arrows_with_ids.id == [1, 2, 3]
    
    def test_arrow_queries(self):
        arrows = DeArrows('data.txt', 'data1.txt', metadata = [{'experiment_name': 'Test Experiment'}, {'experiment_name': 'Test Experiment 2'}], heal_genes = True)
        arrows_experiments = arrows.get_experiment(name = 'Test Experiment')
        assert arrows.get_downregulated().height == 15
        assert arrows.get_upregulated().height == 22
        assert arrows.get_significant_genes().height == 37
        assert len(arrows_experiments.experiment_metadata) == 2
    
    def test_arrow_add(self):
        arrow = DeArrow('data.txt', metadata = {'experiment_name': 'Test Experiment'}, heal_genes = True)
        arrows = arrow.add_experiment('data1.txt', metadata = {'experiment_name': 'Test Experiment 2'}, experiment_id = 3)
        arrows2 = arrows.add_experiment('data.txt', metadata = {'experiment_name': 'Test Experiment 3'}, experiment_id = 4)
        arrows3 = arrows2.add_experiment(arrows, metadata = {'experiment_name': 'Test Experiment 4'}, experiment_id = [7, 8])
        arrows4 = arrows3.add_experiment(arrow, experiment_id = 10)
        assert arrows._table.height > 0
        assert arrows.get_gene('PTEN').height > 0 and arrows.get_gene(ensembl_id = 'ENSMUSG00000000017').height > 0
        assert arrows.id == [1, 3]
        assert arrows2.id == [1, 3, 4]
        assert arrows3.id == [1, 3, 4, 7, 8]
        assert arrows4.id == [1, 3, 4, 7, 8, 10]
    
    def test_arrows_set_id(self):
        arrows = DeArrows('data.txt', 'data1.txt', metadata = [{'experiment_name': 'Test Experiment'}, {'experiment_name': 'Test Experiment 2'}])
        arrows1 = arrows.set_id([4, 5])
        arrows2 = arrows.set_id({2:7})
        assert arrows1.id == [4, 5]
        assert arrows2.id == [1, 7]
    
    def test_arrow_insert(self):
        arrow = DeArrow('data.txt', metadata = {'experiment_name': 'Test Experiment'}, heal_genes = True)
        arrows = arrow.add_experiment('data1.txt', metadata = {'experiment_name': 'Test Experiment 2'}, experiment_id = 3)
        with DeQuackling() as db:
            arrows.insert('SQL.duckdb')
            arrow.insert('SQL.duckdb')
            assert db.get_gene().height > 0
            assert db.get_upregulated().height > 0
    
    def test_write_parquet(self):
        arrow = DeArrow('data.txt', metadata = {'experiment_name': 'Test Experiment'}, heal_genes = True)
        arrows = arrow.add_experiment('data1.txt', metadata = {'experiment_name': 'Test Experiment 2'}, experiment_id = 3)
        arrows.write_parquet(output_path = 'genes.parquet', experiment_id = 1, compression = 'zstd', compression_level = 3)
        assert os.path.exists('genes.parquet')
        assert os.path.exists('genes_metadata.json')
        arrow_parquet = DeArrow('genes.parquet', metadata = {'experiment_name': 'Test Experiment 3'})
        assert arrow_parquet.height > 0
        
    def test_arrow_cleanup(self):
        os.remove('SQL.duckdb')
        os.remove('genes.parquet')
        os.remove('genes_metadata.json')
        
    

