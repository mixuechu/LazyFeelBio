import os
import zipfile
import shutil
import argparse
import logging
import time
from pathlib import Path

# 设置日志记录
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')


def unzip_files(zip_file_paths, temp_dir):
    """
    Unzip all the zip files to a temporary directory.

    Args:
    zip_file_paths (list of str): List of zip file paths.
    temp_dir (str): Temporary directory to store the unzipped files.

    Returns:
    list of str: List of paths to the unzipped files.
    """
    unzipped_files = []
    for zip_file_path in zip_file_paths:
        logging.info(f'Starting to unzip {zip_file_path}')
        with zipfile.ZipFile(zip_file_path, 'r') as zip_ref:
            zip_ref.extractall(temp_dir)
            unzipped_files.extend([os.path.join(temp_dir, f) for f in zip_ref.namelist()])
            logging.info(f'Finished unzipping {zip_file_path}')
        time.sleep(2)  # 模拟延迟
    return unzipped_files


def zip_files(files, output_zip_path):
    """
    Zip the given files into a single zip file.

    Args:
    files (list of str): List of file paths to be zipped.
    output_zip_path (str): Output zip file path.
    """
    logging.info(f'Starting to zip files into {output_zip_path}')
    with zipfile.ZipFile(output_zip_path, 'w') as zipf:
        for file in files:
            arcname = os.path.relpath(file, start=os.path.dirname(files[0]))
            zipf.write(file, arcname)
            logging.info(f'Added {file} to the zip archive')
            time.sleep(1)  # 模拟延迟
    logging.info(f'Finished zipping files into {output_zip_path}')


def main(zip_file_paths, output_zip_path):
    temp_dir = Path('temp_unzipped')
    temp_dir.mkdir(exist_ok=True)

    try:
        # Unzip all files to the temporary directory
        logging.info('Starting the unzip process')
        unzipped_files = unzip_files(zip_file_paths, temp_dir)
        logging.info('Unzip process completed')

        # Create a single zip file containing all unzipped files
        logging.info('Starting the zip process')
        zip_files(unzipped_files, output_zip_path)
        logging.info('Zip process completed')
    finally:
        # Clean up the temporary directory
        logging.info('Cleaning up the temporary directory')
        shutil.rmtree(temp_dir)
        logging.info('Cleanup completed')


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Unzip multiple zip files and recompress all contents into a single zip file.")
    parser.add_argument('zip_files', metavar='ZIP_FILE', type=str, nargs='+', help='Paths to the input zip files.')
    parser.add_argument('-o', '--output', required=True, type=str, help='Path to the output zip file.')

    args = parser.parse_args()
    main(args.zip_files, args.output)