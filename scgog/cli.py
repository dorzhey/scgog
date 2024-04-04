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


def main():
    print("QC, normalization, and Annotation merge")
    data_path = windows_path_parser(input("Please enter the absolute path to your 10x data file in .h5 format: ").strip('"'))
    ann_path = windows_path_parser(input("Please enter the absolute path to your ann data file in .tsv format: ").strip('"'))
    mdata = preprocess_omics_data(data_path, ann_path)
    print("Data preprocessing completed successfully")
    print("Visualization")
    saved = yes_no_answer("Would you like compilated figure to be saved? If not it will be shown")
    file_path = visualize_umap(mdata, saved)
    print("Data visualization completed successfully")
    if saved:
        print("Figure saved at:")
        print(file_path)

    batch_size = 32
    shuffle = yes_no_answer("Would you like Data Loader to shuffle the entries?")
    benchmark_models(mdata, batch_size, shuffle)
    print("ML benchmark completed successfully")

