#!/usr/bin/env python3
from zipfile import ZipFile
from textwrap import dedent
import os
import json
import click
import jsonlines
import re
import sys

class DatasetManager:
    """
    Class to manage and extract data from NCBI genome dataset ZIP files.
    """
    def __init__(self, file: str) -> None:
        self.file: str = file
        self.root_path : str = 'ncbi_dataset/data/'
        self.directory : str = os.path.dirname(self.file)
        self.filename : str = os.path.basename(self.file)
        self.catalog : dict = self.get_catalog()
        self.assembly_report: list = self.get_assembly_report()
        self.genome_path : str = self.get_genome_filepath()
        self.accession : str = self.get_accession()
        self.corrected_accession = self.accession.replace('.', '_')  
        self.genome_extension: str = '.fna'
        self.species: str = self.get_species_name()
        self.genus: str = self.species.split(' ')[0]
        self.strain: str = self.get_strain()
        self.tax_id: str = self.assembly_report[0]['organism']['taxId']
    
    def unzip(self)-> None:
        """Extracts the entire dataset content."""
        with ZipFile(self.file, 'r') as file: 
            file.extractall(path = self.directory)
    
    def get_species_name(self):
        """Retrieves and normalizes the organism name from the assembly report."""
        name = self.assembly_report[0]['organism']['organismName']
        return self.normalize_string(name)

    def normalize_string(self, string: str)-> str:
        """Removes special characters and cleans the species name."""
        result = re.sub(r'[^a-zA-Z0-9-_ ]', '', string)
        return self.remove_strain_from_species(result)
    
    def remove_strain_from_species(self, species: str)-> str:
        """Assumes 'genus species strain' format and keeps only 'genus species'."""
        temp = species.strip().split(' ')
        if len(temp) >= 3:
            return ' '.join(temp[0:2])
        else:
            return species
    
    def unzip_file(self, file: str, output = None, filename = None)-> None:
        """Extracts a single file from the ZIP archive. Robust against null paths."""
        if file is None:
            return
            
        output_path = output if output != None else self.directory
        filename = os.path.basename(file) if filename == None else filename
        
        if not os.path.exists(output_path):
            os.makedirs(output_path)
            
        with ZipFile(self.file) as z:
            try:
                with open(os.path.join(output_path, filename), 'wb') as f:
                    f.write(z.read(file))
            except (KeyError, Exception) as e:
                # Extraction errors are sent to stderr to avoid breaking the kSNP index on stdout
                print(f"Warning: Could not extract {file} from {self.filename}: {e}", file=sys.stderr)
    
    def load_zipped_file(self, file_inside_zip: str)-> list:
        """Returns the content of a file located inside the ZIP archive."""
        with ZipFile(self.file) as myzip:
            with myzip.open(file_inside_zip) as myfile:
                result = myfile.readlines()
                return list(result)
    
    def load_zipped_jsonlines(self, file_inside_zip)-> list:
        """Returns a list of dictionaries from an internal JSONL file."""
        bytes_content = self.load_zipped_file(file_inside_zip)
        result = []
        for line in jsonlines.Reader(bytes_content):
            result.append(line)
        return result

    def get_filepath(self, file_type: str)-> str:
        """Retrieves the internal path for a file type based on the catalog."""
        legend = {
            'CDS':'CDS_NUCLEOTIDE_FASTA',
            'GENOME':'GENOMIC_NUCLEOTIDE_FASTA',
            'GFF':'GFF3',
            'REPORT':'SEQUENCE_REPORT'
        }
        try:
            files = self.catalog['assemblies'][1]['files']
            for f in files:
                if f['fileType'] == legend[file_type]:
                    path = f['filePath']
                    return f'{self.root_path}{path}' 
        except (KeyError, IndexError):
            pass
        return None

    def get_catalog(self)-> dict:
        """Parses the dataset_catalog.json file."""
        bytes_content = self.load_zipped_file(self.root_path + 'dataset_catalog.json')
        result = json.loads(''.join([i.decode('utf-8') for i in bytes_content]))
        return result
    
    def get_assembly_report(self)-> list:
        """Parses the assembly_data_report.jsonl file."""
        path = self.root_path + 'assembly_data_report.jsonl'
        return self.load_zipped_jsonlines(path)
        
    def get_genome_filepath(self)-> str:
        return self.get_filepath(file_type='GENOME')
    
    def get_gff_filepath(self)-> str:
        return self.get_filepath(file_type='GFF')
    
    def get_CDS_filepath(self)->str:
        return self.get_filepath(file_type='CDS')

    def get_accession(self)-> str:
        return self.catalog['assemblies'][1]['accession']
    
    def get_strain(self)->str:
        """Extracts strain name if available, otherwise returns 'NA'."""
        try:
            result = self.assembly_report[0]['organism']['infraspecificNames']['strain']
        except KeyError:
            result = 'NA'
        return result

    def extract_assembly_report(self, output_path = None)-> None:
        path_to_report = self.root_path + 'assembly_data_report.jsonl'
        self.unzip_file(path_to_report, output_path)
    
    def extract_genome(self, output_path = None)-> None:
        filename = f'GENOME_{self.corrected_accession}{self.genome_extension}'
        self.unzip_file(self.genome_path, output_path, filename)
    
    def extract_CDS(self, output_path = None)-> None:
        filename = f'CDS_{self.corrected_accession}{self.genome_extension}'
        cds_path = self.get_CDS_filepath()
        self.unzip_file(cds_path, output_path, filename)
    
    def extract_gff(self, output_path = None)-> None:
        path_to_gff = self.get_gff_filepath()
        if path_to_gff:
            filename = 'genomic.gff'
            self.unzip_file(path_to_gff, output_path, filename)

    def build_dataset_dirname(self):
        return self.filename

    def save_organism_file(self, filename)-> None:
        """Saves metadata to a TSV file for record keeping."""
        content = dedent(f'''\
            genus\t{self.genus}
            species\t{self.species}
            strain\t{self.strain}
            tax-id\t{self.tax_id}
            accession\t{self.accession}
            genome\tGENOME_{self.corrected_accession}{self.genome_extension}
        ''')
        with open(filename, mode='w') as file:
            file.write(content)

    def build_dataset(self, output_path = None)-> None:
        """
        Orchestrates the uncompressed dataset construction.
        Outputs a tab-separated report for kSNP to stdout.
        """
        output = output_path if output_path != None else self.directory
        dirname = self.build_dataset_dirname()
        dataset_path = os.path.join(output, os.path.splitext(dirname)[0])
        
        # 1. Extract Genome (High Priority)
        try:
            self.extract_genome(output_path = dataset_path)
        except Exception as e:
            print(f"Fatal error extracting genome {self.filename}: {e}", file=sys.stderr)
            return

        # 2. Secondary extractions (Non-blocking)
        try: self.extract_CDS(output_path = dataset_path)
        except: pass
        
        try: self.extract_assembly_report(output_path = dataset_path)
        except: pass
        
        try: self.extract_gff(output_path = dataset_path)
        except: pass

        # 3. Save Metadata
        try:
            organism_filepath = os.path.join(dataset_path, 'organism_report.tsv')
            self.save_organism_file(organism_filepath)
        except: pass

        # 4. Final report for kSNP index (Stdout)
        genome_type = 'GENOME'
        genome_path = os.path.normpath(
            f'{dataset_path}/{genome_type}_{self.corrected_accession}{self.genome_extension}'
        )
        report = [self.genus, self.species, self.strain, str(self.tax_id), genome_path]
        click.echo('\t'.join(report))

# CLI Configuration
@click.group()
def cli():
    pass

@cli.command()
@click.argument('ncbi_dataset')
@click.option('--output', default = None, help = 'Output path for FASTA files')
def extract_fasta(ncbi_dataset: str, output: str)-> None:
    """Extracts only the genomic FASTA file."""
    dataset = DatasetManager(ncbi_dataset)
    dataset.extract_genome(output_path = output)

@cli.command()
@click.argument('ncbi_dataset')
@click.option('--output', default = None, help = 'Directory to build the uncompressed dataset')
def build_dataset(ncbi_dataset: str, output: str)-> None:
    """Extracts all relevant files and prints the metadata report."""
    dataset = DatasetManager(ncbi_dataset)
    dataset.build_dataset(output_path = output)

if __name__ == '__main__':
    cli()