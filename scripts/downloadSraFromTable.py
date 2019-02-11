#!/usr/bin/env python
import shutil, os, argparse, sys, stat, urllib

from urlgrabber.grabber import URLGrabber
import requests
import csv
from joblib import Parallel, delayed
import traceback

def get_file_if_size_diff(url, d):
    fn = url.split('/')[-1]
    out_fnp = os.path.join(d, fn)
    g = URLGrabber(reget = "simple")
    locFnp = g.urlgrab(url, out_fnp)
    return locFnp

def downloadFileAttempt(url, directoryName):
    try:
        return get_file_if_size_diff(url, directoryName)
    except Exception as err:
        print ("Failed  to download " + str(url) + " to " + str(directoryName) + ", mess: " + str(err))
        traceback.print_exc()

def parse_args_downloadFiles():
    parser = argparse.ArgumentParser()
    parser.add_argument( type=str, required = True, help = "SRA url file created by getSraRunsFromAccIds.py")
    parser.add_argument('-o','--outDirectory', type=str, default = "./", help = "The directory in which to download the files")
    parser.add_argument('-n','--ncpus', type=int, default = 1, help = "Number of cpus to use")
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args_downloadFiles()
    with open(args.sraUrlFnp, 'rb') as urlFile:
        reader =  csv.DictReader(urlFile, delimiter='\t')
        urls = []
        for row in reader:
            urls.append(row.get("url"))
        downloadedFiles = Parallel(n_jobs = args.ncpus)(delayed(downloadFileAttempt)(url, "./") for url in urls )
        print(downloadedFiles)
