from data_preprocessing import *
from visualization import *
import pathlib

def windows_path_parser(path):
    return pathlib.PureWindowsPath(path).as_posix()

def yes_no_answer(input_msg, end_msg, func, param):
    while True:

        choice = input(input_msg+" [y/n] ")
        if choice == 'y':
            result =  func(param)
            print(end_msg)
            return result
        elif choice == 'n':
            return param
        else:        
            print("incorrect input, try again")


def main():
    data_path = windows_path_parser(input("Please enter the aboslute path to your scRNA-seq data file in .h5 format: "))
    mdata = preprocess_omics_data(data_path)
    print("Data input completed successfully")

    mdata = yes_no_answer("Do you want to QC the data?","Data quality control completed successfully", quality_control, mdata)
    file_path = yes_no_answer("Do you want to visualize the data?","Here are the plots", visualize_umap, mdata)
    print(f"UMAP plot saved to: {file_path}")
    mdata = yes_no_answer("Do you want to intersect the data?","Merge completed", merge_data, mdata)
    

main()
    