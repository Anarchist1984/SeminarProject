from src.chembl import ChemBLDataProcessor
uniprot_id = "P14780"
processor = ChemBLDataProcessor(uniprot_id)
processed_data = processor.request_similar_compounds()
processed_data.head()