#!env python2
from __future__ import print_function
import argparse
import sys
import pandas as pd
import time
import os
import glob
import astropy.io.fits as pyfits
import h5py
import random
import warnings
import traceback
import numpy as np

from tilemaker import find_tiles, make_cuts
from getfile import download_file
from shutil import copyfile
from multiprocessing import Process, Manager

warnings.simplefilter('ignore')
MAX_TILE_ATTEMPTS = 3

verbose = [False]
delete_fits = [True]

symbols = "!@#$%^&*()-=_+[]{}\|,./<>?"

tiles = pd.read_csv(
    os.path.dirname(os.path.realpath(__file__)) + "/y3a1tiles.csv",
    index_col="TILENAME")


def progbar(current, to, width=40, show=True, message=None):
    try:
        percent = float(current) / float(to)
    except:
        percent = 0
    length = int(width * percent)
    if show:
        count = " (%d/%d)    " % (current, to)
    else:
        count = ""
    if message:
        count += message
    sys.stdout.write(("\r[" + ("#" * length) + " " * (width - length) + "] %0d"
                      % (percent * 100)) + "%" + count)
    sys.stdout.flush()


class StampWorker(object):
    def __init__(self, rank, name, catalog, datastore, datamap, dimension):
        self.datamap = datamap
        self.rank = rank
        self.name = name
        self.catalog = catalog
        self.datastore = datastore
        self.done = datamap['done_tiles']
        self.tile_list = list(catalog.TILENAME.unique())
        self.problem_tiles = {}
        self.complete_tiles = []
        self.dimension = dimension

    def printlog(self, logstring):
        if verbose[0]:
            print(logstring)
        self.log_to_file(logstring)

    def log_to_file(self, logstring):
        with open("logfile_" + str(self.rank) + ".log", "a") as f:
            f.write(logstring + "\n")

    def get_next_tile(self):
        return self.tile_list.pop()

    def run(self, delay):
        datamap = self.datamap
        time.sleep(delay)
        self.printlog("Starting " + self.name)
        NN = len(self.tile_list)
        while len(self.tile_list) > 0:
            try:
                tile_to_process = self.get_next_tile()

                if self.rank == 0:
                    progbar(len(self.done), NN)

                if tile_to_process is None:
                    break
                self.printlog(self.name + " working on tile " +
                              tile_to_process)
                filenames = get_tile_files(tile_to_process)
                if filenames is None:
                    self.printlog("No files to process for tile %s. Abort!" %
                                  tile_to_process)
                    continue
                tile_sources = self.catalog[self.catalog['TILENAME'] ==
                                            tile_to_process]
                ds = self.datastore["/stamps/" + tile_to_process + "/data"]
                masks = self.datastore["/stamps/" + tile_to_process + "/masks"]
                hd = self.datastore["/stamps/" + tile_to_process + "/header"]
                catmeta = self.datastore['/stamps/' +
                                         tile_to_process + "/catalog"]
                for f in filenames:
                    self.printlog("%d Making cuts to %s, %s" %
                                  (self.rank, str(f[0]), str(f[1])))

                    with warnings.catch_warnings():
                        warnings.simplefilter("ignore")
                        fits = pyfits.open(f[0])

                    results = make_cuts(
                        tile_sources,
                        fits,
                        f[1],
                        self.dimension,
                        tofile=False,
                        data=ds,
                        headers=hd,
                        catmeta=catmeta,
                        masks=masks,
                        logfile="logfile_" + str(self.rank) + ".log")

                    self.printlog("Made cuts successfully.")
                    datamap['bad_objects'].extend(results['bad_objects'])

                self.done.append(tile_to_process)
                fits.close()
                cleanup(tile_to_process)
                self.printlog("Done with tile " + tile_to_process)
            except:
                self.printlog("Error in main loop of worker %d. Tile: %s" %
                      (self.rank, tile_to_process))
                self.printlog("Unexpected error: " + str(sys.exc_info()[0]))
                traceback.print_exc(file=sys.stdout)

                # Try tile again
                attempts = self.problem_tiles.get(tile_to_process, 1)
                if attempts < MAX_TILE_ATTEMPTS:
                    self.tile_list.append(tile_to_process)
                    self.problem_tiles[tile_to_process] = attempts + 1
                    self.printlog("Will try again. Attempts: %d" % attempts)
                else:
                    datamap['failed_tiles'].append(tile_to_process)

        if self.rank == 0:
            progbar(NN, NN)


def main(argv):
    """Usage:

    catalog_to_stamps.py -p <procs> """
    parser = argparse.ArgumentParser()

    # Options
    parser.add_argument(
        '--processes', type=int, default=1, help="Number of cores to use.")
    parser.add_argument(
        '--dimension',
        type=int,
        default=100,
        help="Size of stamps in pixels, default=100")
    parser.add_argument(
        '--flatten',
        help="Don't store stamps by tile in HDF5 output.",
        action="store_true")
    parser.add_argument(
        '--no-cleanup',
        help="Don't remove fits files after processing.",
        action="store_true")

    parser.add_argument(
        '--verbose',
        help="More logging info to console.",
        action="store_true")
    # Positional arguments
    parser.add_argument(
        "input_catalog",
        help="Catalog containing RA, DEC of objects to extract")
    parser.add_argument(
        "tiles_per_batch",
        type=int,
        help="Number of survey tiles to include in each output data file.")
    parser.add_argument("batches", type=int, help="Number of batches to run.")
    parser.add_argument("datastore_prefix", help="Name of output data file.")
    args = parser.parse_args()

    dimension = args.dimension

    if dimension > 100:
        print("Warning: Making stamps of size %dx%d." % (dimension, dimension))
    if args.verbose:
        verbose[0] = True
    if args.no_cleanup:
        delete_fits[0] = False

    catalog_file = args.input_catalog
    batch_size = args.tiles_per_batch
    batches = args.batches
    dstore_prefix = args.datastore_prefix

    procs = args.processes
    start_index = 1
    while True:
        if os.path.exists(dstore_prefix + "_" + str(start_index) + ".hdf5"):
            start_index += 1
        else:
            break

    copyfile(catalog_file, catalog_file + ".bak")
    catalog = pd.read_csv(
        catalog_file, quotechar='"', index_col='COADD_OBJECT_ID')
    catalog = catalog[catalog.TILENAME != "NONE"]

    catalog_toproc = catalog[catalog.STATUS == 'new']
    tilegroups = catalog_toproc.groupby(by="TILENAME")
    tilegroups = tilegroups.size().sort_values(ascending=False)
    tilenames = tilegroups.index.values

    # Start the processes
    manager = Manager()
    results_dict = {}
    b = 0
    while b < batches:
        print("\nBatch %d/%d." % (1 + (b/procs), batches))
        proc_handles = []
        for p in range(procs):
            if b >= batches:
                break
            batch_dict = manager.dict()
            dstore_name = "%s_%0.3d.hdf5" % (dstore_prefix, start_index)
            results_dict[dstore_name] = batch_dict
            start = b * batch_size
            # start = ((b * p) + p) * batch_size
            end = start + batch_size
            tilebatch = catalog_toproc[catalog_toproc.TILENAME.isin(tilenames[
                start:end])]

            try:
                proc = Process(
                    target=main_batch,
                    args=(tilebatch, dstore_name, p, args.flatten, dimension,
                          batch_dict))
                proc.start()
                proc_handles.append(proc)
                b += 1
                start_index += 1
            except StopIteration:
                break
        while proc_handles:
            proc = proc_handles.pop()
            proc.join()

    for batch in results_dict.keys():
        batch_results = results_dict[batch]
        to_update = batch_results['done_tiles']
        failed_tiles = batch_results['failed_tiles']
        bad_objects = batch_results['bad_objects']

        catalog.loc[catalog.TILENAME.isin(to_update), "STATUS"] = 'done'
        catalog.loc[catalog.TILENAME.isin(failed_tiles), "STATUS"] = 'failed'
        catalog.loc[bad_objects, "STATUS"] = 'failed'
    catalog.to_csv(catalog_file)


def main_batch(catalog, dstore, rank, flatten, dimension, results_dict):
    if verbose[0]:
        print("\nMaking stamps from %d tiles into datastore %s on core %d" %
          (len(catalog.TILENAME.unique()), dstore, rank+1))
    if not os.path.exists(dstore):
        initialise_datastore(dstore, dimension)

    datastore = h5py.File(dstore, 'r+')
    results_dict['failed_tiles'] = []
    results_dict['done_tiles'] = []
    results_dict['bad_objects'] = []

    N = len(catalog.index)
    seconds = random.randint(0, 2)
    time.sleep(seconds)

    if flatten:
        datastore.create_dataset(
            "/stamps/data", (N, dimension, dimension, 5),
            dtype='f4')
        datastore.create_dataset(
            "/stamps/masks", (N, 5), dtype=np.int32)
        datastore.create_dataset(
            "/stamps/header", (N, 5), dtype="S9000")
        datastore.create_dataset(
            "/stamps/catalog", (N, ), dtype='S30')

    tilegroups = catalog.groupby(by="TILENAME")
    tilegroups = tilegroups.size().sort_values(ascending=False)
    tilenames = tilegroups.index.values


    for i in range(len(tilenames)):
        tilename = tilenames[i]
        groupsize = tilegroups.iloc[i]
        if not flatten:
            datastore.create_dataset(
                "/stamps/%s/data" % (tilename),
                (groupsize, dimension, dimension, 5),
                dtype='f4')
            datastore.create_dataset(
                "/stamps/%s/masks" % (tilename), (groupsize, 5),
                dtype=np.int32)
            datastore.create_dataset(
                "/stamps/%s/header" % (tilename), (groupsize, 5),
                dtype="S9000")
            datastore.create_dataset(
                "/stamps/%s/catalog" % (tilename), (groupsize, ), dtype='S30')

    datamap = {"bad_objects": [], "failed_tiles": [], "done_tiles": []}
    worker = StampWorker(rank, "Worker" + str(rank), catalog, datastore,
                         datamap, dimension)
    worker.run(random.randint(1, 10))
    datastore.close()

    results_dict['failed_tiles'] += datamap['failed_tiles']
    results_dict['bad_objects'] += datamap['bad_objects']
    results_dict['done_tiles'] += datamap['done_tiles']


def initialise_datastore(datastore, dimension):
    f = h5py.File(datastore, 'w')
    grp = f.create_group("stamps")
    grp.attrs['description'] = "DES y3a1coadd cutouts, dimensions=%dx%d" % (
        dimension, dimension)


def get_tile_files(tile):
    try:
        tilefile = tiles.loc[tile, 'FILENAME']
        tilepath = tiles.loc[tile, 'PATH']
    except ValueError:
        print("No file for tile " + tile + ". Abort!")
        return None

    if type(tilefile) != str:  # Use first one from list
        tilefile = tilefile[0]
        tilepath = tilepath[0]

    filenames = []
    for band in "grizY":
        filetoget = tilefile.replace("_r.fits", "_" + band + ".fits.fz")
        filenames.append((filetoget, band))
        if not os.path.isfile(filetoget):
            download_file(tilepath + "/" + filetoget)
    return filenames


def cleanup(tile):
    fitsfiles = glob.glob(tile + "*.fz")
    # print("Deleting " + str(fitsfiles))
    for ff in fitsfiles:
        if delete_fits[0]:
            os.remove(ff)


def to_tiles(catalog, output):
    find_tiles(catalog)
    catalog['STATUS'] = 'new'
    if 'COADD_OBJECT_ID' not in catalog.columns:
        catalog['COADD_OBJECT_ID'] = catalog['name']
    catalog.set_index(['COADD_OBJECT_ID'])

    catalog.to_csv(output)

    groups = catalog.groupby(by="TILENAME")
    print("Sources are on %d tiles with an average of %.2f per tile." %
          (groups.ngroups, groups.size().mean()))
    print("Top tiles by number of sources: ")
    print(groups.size().sort_values(ascending=False).head(10))

    # tiles = groups.name


if __name__ == "__main__":
    start = time.time()
    main(sys.argv)
    print("\nDone in %.2f minutes." % ((time.time() - start) / 60.0))
