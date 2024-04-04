from data_preprocessing import *
from visualization import *
import pathlib
from benchmarking_ml_models import *

def windows_path_parser(path):
    return pathlib.PureWindowsPath(path).as_posix()

def yes_no_answer(input_msg):
    for _ in range(3):

        choice = input(input_msg+" [y/n] ")
        if choice == 'y':
            return True
        elif choice == 'n':
            return False
        else:        
            print("incorrect input, try again")


def preprocess():
    print("QC, normalizatio, and Annotation merge")
    data_path = windows_path_parser(input("Please enter the absolute path to your 10x data file in .h5 format: ").strip('"'))
    ann_path = windows_path_parser(input("Please enter the absolute path to your ann data file in .tsv format: ").strip('"'))
    save_path = preprocess_omics_data(data_path, ann_path)
    print("Data preprocessing completed successfully and MuData object saved at:")
    print(save_path)

def visualize():
    mdata_path = windows_path_parser(input("Please enter the absolute path to your MuData file in .h5mu format: ").strip('"'))
    saved = yes_no_answer("Would you like compilated figure to be saved? If not it will be shown")
    file_path = visualize_umap(mdata_path, saved)
    print("Data visualization completed successfully")
    if saved:
        print("Figure saved at:")
        print(file_path)

def ml_benchmarking():
    mdata_path = windows_path_parser(input("Please enter the absolute path to your MuData file in .h5mu format: ").strip('"'))
    label = yes_no_answer("Would you like to provide name of target variable? If not 'gene' will be used")
    if label:
        label = input("Please type name of target as it appears in Annotation file")
    else:
        label = "gene"
    batch_size = yes_no_answer("Would you like to provide name of batch size? If not 32 will be used")
    if batch_size:
        batch_size = input("Please type name of target as it appears in Annotation file")
    else:
        batch_size = 32
    shuffle = yes_no_answer("Would you like Data Loader to shuffle the entries?")
    file_path = benchmark_models(mdata_path,label, batch_size, shuffle)
    print("ML benchmark completed successfully")

