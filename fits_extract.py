#!env python
import os
import sys
import coltools as ct
import h5py as h5
import astropy.io.fits as pyfits


def fal(filename, skip_comments=True):
    with open(filename, "r") as f:
        lines = f.readlines()
    if skip_comments:
        return [l.strip() for l in lines if not l.startswith("#")]
    else:
        return [l.strip() for l in lines]


def progbar(current, to, width=40, show=True, message=None, stderr=False):
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


def extract_objects(datastore, objids, outdir):
    if objids is None:
        extract_all(datastore, outdir)
        return
    else:
        extract_some(objids, datastore, outdir)


def extract_some(ids, datastore, outdir):
    todo = len(ids)
    done = 0
    ct.progbar(done, todo)
    with h5.File(datastore, "r") as f:
        ds = f['/stamps/']
        for tile in ds:
            data = ds[tile + "/data"]
            headers = ds[tile + "/header"]
            objids = ds[tile + "/catalog"]
            N = data.shape[0]
            for i in range(N):
                objid = objids[i].strip()
                if objid in ids:
                    done += 1
                    ct.progbar(done, todo)
                    heads = headers[i]
                    for b in range(5):
                        band_name = "grizY" [b]
                        d = data[i, :, :, b]
                        head = pyfits.Header.fromstring(heads[b].strip())
                        filename = "%s/%s_%s.fits" % (outdir, objid, band_name)
                        outfits = pyfits.PrimaryHDU(data=d, header=head)
                        outfits.writeto(filename, overwrite=True)


def make_rgb(datastore, outdir):
    pass


def extract_all(datastore, outdir):
    with h5.File(datastore, "r") as f:
        ds = f['/stamps/']
        for tile in ds:
            data = ds[tile + "/data"]
            headers = ds[tile + "/header"]
            objids = ds[tile + "/catalog"]
            N = data.shape[0]
            for i in range(N):
                objid = objids[i].strip()
                heads = headers[i]
                for b in range(5):
                    band_name = "grizY" [b]
                    d = data[i, :, :, b]
                    head = pyfits.Header.fromstring(heads[b].strip())
                    filename = "%s/%s_%s.fits" % (outdir, objid, band_name)
                    outfits = pyfits.PrimaryHDU(data=d, header=head)
                    outfits.writeto(filename, overwrite=True)


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("image_extract <datastore> <dest_dir> [objids]")
        sys.exit(0)

    datastore = sys.argv[1]
    outdir = sys.argv[2]

    inputs = None
    if len(sys.argv) == 4:
        inputs = sys.argv[3]

    if not os.path.exists(outdir):
        print("Output directory doesn't exist.")
        sys.exit(0)


    if inputs is not None and os.path.exists(inputs):
        objids = ct.fal(inputs)
    elif inputs is not None:
        objids = inputs.split(",")
    else:
        objids = None

    print(objids)
    extract_objects(datastore, objids, outdir)
