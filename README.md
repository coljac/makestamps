# Making DES stamps

This code will extract a catalog of postage stamps from the DES Y3A1 tiles. It sorts your catalog of targets by DES tile then iterates through all the tiles - all 10,000 - makes the cutouts and stores them in HDF5 files for portability.

## Prerequisites

The requirements can be installed with:

    pip install -r requirements.txt

Basically it needs:
- Python 2.7 (May work with 3 - will test)
- Astropy
- Numpy
- requests to grab the data

(You'll want Cython too.)

## Prepare the catalog

The stamp extractor expects a comma-delimited file with the following headers:

- RA
- DEC
- COADD_OBJECT_ID: This doesn't have to be a real DES id, but this will be the label associated with the fits file.
- TILENAME: The des tile on which the object is to be found.
- STATUS: initial value must be "new"; the program will record the stamps made ("done") or errors here, so repeated invocations of the script won't download the same stamps again by default.

If you don't have the tile data already, rather than going back to the Oracle DB, run the following to add this information to the catalog as a pre-processing step:

    python tilemaker.py <input catalog> <output file>

For example:

    python tilemaker.py my_catalog.csv catalog_with_tiles.csv

For this to work, you'll need a file containing the information on the DES tiles (y3a1tiles.csv) which currently isn't in this repo.

## Run the stamp-maker

The stamp-maker will process all the stamps on the supplied number of tiles and put the stamps into a HDF5 file for local storage. When running the script, the key parameters are:

- Tiles per batch: How many tiles worth of stamps to put in an output files.
- Number of batches: How many output files to create.
- Number of cores to use: If this is greater than one, it will process several batches in parallel.

    python catalog_to_stamps.py --processes <num of threads> --dimension <cutout size> [catalog file] [tiles_to_process] [tiles_per_output] [file prefix]

For example:

    python catalog_to_stamps.py --processes 4 --dimension 128 mycat.csv 100 20 mystamps

Will take the catalog file `mycat.csv` and get the stamps from the top DES tiles, putting 100 tiles worth of stamps in each of 20 hdf5 outputs, named mystamps_001.hdf5 to mystamps_020.hdf5. It will run in parallel on 4 cores.

For more info:

    python catalog_to_stamps.py --help

## Getting individual FITS files

To extract individual FITS stamp files for easier use, use the h5tofits utility:

    python fits_extract.py <HDF5 file> <output dir>

or

    python fits_extract.py <HDF5 file> <output dir> <list_of_objects> 

The list of objects can be either:
- A text file with one object id per line
- A comma-delimited string: object1,[object2,object3,...objectN]

For example:

    python <list_of_objects> mycat_001.hdf5 outputs/ 1234567,

## Structure of output file

The HDF5 files output by the stamp maker contain the following datasets, one per tile processed:

/stamps/<tile>/data: FITS data; dimensions N x dim x dim x 5, where N = stamps in this tile, dim = image dimensions, 5 = number of bands.
/stamps/<tile>/masks: The number of masked pixels according to the mask HDU in the source FITS file.
/stamps/<tile>/header: FITS headers, N x 5
/stamps/<tile>/catalog: Object ids corresponding to the data, dimensions N x 1

## Sample data

To test the stamp maker, a sample catalog is included.

    python catalog_to_stamps.py --processes 1 --dimension 100 sample.csv 1 1 sample

This should produce the output file sample_001.hdf5.

To extract fits files individually:

    mkdir fits
    python fits_extract.py sample_001.hdf5 fits

There should be 245 fits stamps in the fits/ directory.

The following would also work:

    python fits_extract.py sample_001.hdf4 fits 140013384,140013438

And creates 10 files, five bands for each of two objects.
