import numpy as np
import pandas as pd
import h5py
import astropy
import time
import os
import sys
import astropy.io.fits as pyfits
from astropy.coordinates import SkyCoord
from astropy.wcs import WCS
from astropy.nddata import Cutout2D
from astropy import units as u

tiles = pd.read_csv(
    os.path.dirname(os.path.realpath(__file__)) + "/y3a1tiles.csv")


def pb(current, to, width=40, show=True, message=None, stderr=False):
    percent = float(current) / float(to)
    length = int(width * percent)
    if show:
        count = " (%d/%d)    " % (current, to)
    else:
        count = ""
    if message:
        count += message
    outstream = sys.stderr if stderr else sys.stdout
    outstream.write(("\r[" + ("#" * length) + " " * (width - length) +
                     "] %0d" % (percent * 100)) + "%" + count)
    outstream.flush()


def grab_tile(tilename, bands=None):
    pass


def log_to_file(logfile, logstring):
    if not logfile:
        return
    with open(logfile, "a") as f:
        f.write(logstring + "\n")


def make_cuts(catalog,
              tile,
              band,
              stamp_size,
              tofile=True,
              data=None,
              headers=None,
              catmeta=None,
              masks=None,
              logfile=None,
              results=None):
    results = dict(bad_objects=[])
    """Given a fits data file, turn WCS into pixels and grab data from tile."""
    band_idx = "grizY".index(band)
    todo = len(catalog)
    i = 0
    w = WCS(tile[1].header)
    cutouts = []
    mask_sums = []
    filenames = []
    objids = []
    log_to_file(logfile, "Starting cutouts with tile ")

    for obj in catalog[['RA', 'DEC']].itertuples():
        # pb(i + 1, todo)
        # _, objid, ra, dec = obj
        objid, ra, dec = obj
        x, y = w.all_world2pix(catalog.iloc[i]['RA'], catalog.iloc[i]['DEC'],
                               1)
        mask_sum = 0
        try:
            cutout = Cutout2D(
                tile[1].data, (x, y), (stamp_size, stamp_size), wcs=w)
            if masks is not None:
                mask_cut = Cutout2D(
                    tile[2].data, (x, y), (stamp_size, stamp_size), wcs=w)
                mask_sum = mask_cut.data.sum()
            # print("Mask:", objid, mask_sum)
        except TypeError as e:
            results['bad_objects'].append(objid)
            print("Error with object " + str(objid))
            print(getattr(e, 'message', repr(e)))
            log_to_file(logfile, "Error with object %s" % (str(objid)))
            time.sleep(10)  # In case system overloaded
            continue

        cutouts.append(cutout)
        if masks is not None:
            mask_sums.append(mask_sum)
        filenames.append("stamps/" + str(objid) + "_" + band + ".fits")
        objids.append(objid)
        # write_cut(cutout, "stamps/" + str(objid) +  "_" + band + ".fits", tile)
        i += 1
        if i % 50 == 0:
            log_to_file(logfile,
                        "Done %d cutouts of %d. Still running." % (i, todo))
    log_to_file(logfile, "Done with the cutouts, now to store.")
    # start = time.time()

    for j, cutout in enumerate(cutouts):
        head = tile[1].header.copy()
        head['CRPIX1'] = cutout.wcs.wcs.crpix[0]
        head['CRPIX2'] = cutout.wcs.wcs.crpix[1]

        if tofile:
            write_cut(cutout, filenames[j], tile, head)
        else:
            if cutout.data.shape[0] != stamp_size or cutout.data.shape[1] != stamp_size:
                log_to_file(logfile, "Size mismatch: %d, %d (%d)" % \
                        (cutout.data.shape[0], cutout.data.shape[1], stamp_size))
                zeros = np.zeros((stamp_size, stamp_size))
                zeros[0:cutout.data.shape[0], 0:cutout.data.shape[
                    1]] = cutout.data
                cutout.data = zeros

            data[j, :, :, band_idx] = cutout.data
            if masks is not None:
                masks[j, band_idx] = mask_sums[j]
            headers[j, band_idx] = head.tostring().ljust(9000, ' ')
            catmeta[j] = str(objids[j]).ljust(30, ' ')
            # print(" %d/%d " % (j, todo))
        if j % 50 == 0:
            log_to_file(logfile,
                        "Stored %d cutouts of %d. Not dead yet." % (j, todo))
    log_to_file(logfile, "Complete.")
    # end = time.time()
    # print("Time per write: %.4fs" % ((end-start/len(filenames))))
    return results


def extract_tile(datastore, tile):
    pass


def extract_stamp(datastore, tile, objid, filename):
    write_cut()


def write_cut(cutout, filename, tile, head):
    """Save a cutout to the filesystem as a fits file."""
    outfits = pyfits.PrimaryHDU(data=cutout.data, header=head)
    outfits.writeto(filename, overwrite=True)


def find_tiles_reverse(catalog):
    """Lookup the tile name given an RA, DEC pair."""

    N = len(tiles.index)
    i = 0
    catalog['TILENAME'] = "NONE"
    catalog['STATUS'] = "new"
    for tile in tiles[["URAMIN", "URAMAX", "UDECMIN", "UDECMAX",
                       "TILENAME"]].itertuples():
        pb(i + 1, N)
        idx, ramin, ramax, decmin, decmax, tilename = tile

        if ramin > ramax:
            found = catalog[((catalog.RA > ramin) | (catalog.RA < ramax)) &
                            (catalog.DEC > decmin) & (catalog.DEC < decmax)]
        else:
            found = catalog[(catalog.RA > ramin) & (catalog.RA < ramax) & \
                (catalog.DEC > decmin) & (catalog.DEC < decmax)]
        catalog.loc[found.index, 'TILENAME'] = tilename
        i += 1


def find_tiles(catalog):
    """For a list of ra/dec pairs, find a tile that they are located in."""
    i = 0
    length = len(catalog)
    catalog['TILENAME'] = "---"
    last_tile = None
    for row in catalog[['RA', 'DEC']].itertuples():
        pb(i + 1, length)
        idx, ra, dec = row
        if last_tile is not None and (last_tile.URAMIN < ra) and (last_tile.URAMAX > ra) and\
                (last_tile.UDECMIN < dec) and (last_tile.UDECMAX > dec):
            catalog.loc[idx, 'TILENAME'] = last_tile.TILENAME
            continue

        found = tiles[(tiles.URAMIN < ra) & (tiles.URAMAX > ra) &
                      (tiles.UDECMIN < dec) & (tiles.UDECMAX > dec)]
        if len(found) == 0:
            catalog.loc[idx, 'TILENAME'] = "NONE"
        else:
            catalog.loc[idx, 'TILENAME'] = found.iloc[0].TILENAME
            last_tile = found.iloc[0]
        i += 1


if __name__ == "__main__":
    """Augment an object catalog with DES tile names and a status.

    Usage: tilemaker.py <input cat> <output filename>
    """
    catalog = pd.read_csv(sys.argv[1], index_col="COADD_OBJECT_ID")
    find_tiles_reverse(catalog)
    catalog.to_csv(sys.argv[2])
