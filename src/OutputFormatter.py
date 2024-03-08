import csv
from importlib.resources import files
from typing import Dict, List, Tuple
from click.types import File
from openpyxl import Workbook, load_workbook

from DataClass import Mutation, Sample


def matrix_format(mutations: List[Mutation]) -> Tuple[List[str], List[Dict[str, str]]]:
    header = ['Sample']
    samples = {}
    for mutation in mutations:
        if mutation.found and mutation.name not in header:
            header.append(mutation.name)
        for sample in mutation.samples:
            if sample not in samples:
                samples[sample] = {'Sample': sample}
            samples[sample][mutation.name] = mutation.samples[sample]
    return header, samples.values()


def tabular_output(samples: List[Sample]) -> Tuple[List[str], List[Dict[str, str]]]:
    header = ['Sample', 'Marker mutations', 'Found mutations', 'Effect', 'Subtype', 'Papers']
    data = []
    for sample in samples:
        for marker in sample.markers:
            marker['Sample'] = sample.name
            data.append(marker)
    return header, data


def write_csv(output_file: File, header: List[str], data: List[Dict[str, str]]) -> None:
    writer = csv.DictWriter(output_file, header, delimiter='\t', lineterminator='\n', extrasaction='ignore')
    writer.writeheader()
    writer.writerows(data)    

def get_workbook() -> Workbook:
    return Workbook()

def save_workbook(wb: Workbook, save_path: str) -> None:
    wb.save(save_path)

def write_excel(wb: Workbook, sheet: str, header: List[str], data: List[Dict[str, str]]) -> Workbook:
    wb.create_sheet(sheet)
    ws = wb[sheet]
    for col, value in enumerate(header):
        ws.cell(row=1, column=col+1, value=value)
    for row, values in enumerate(data):
        for col, col_name in enumerate(header):
            ws.cell(row=row+2, column=col+1, value=values.get(col_name, ''))
    return wb
