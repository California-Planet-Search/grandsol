#!/usr/bin/env python

# Import packages.
import argparse
import struct
import numpy as np
from astropy.io import fits

def makehead(liststr):
    """Interpret list of byte-array strings as a FITS header.
    Currently, each byte array must yield a string with 80 characters.
    """
    headstr = ''
    for str in liststr:
        card = str.decode('utf-8')
        assert len(card) == 80
        headstr += card
    return fits.Header().fromstring(headstr)

class Rechead:
    """Contents of a record header from a 'dsk' file.
    """

    def __init__(self, recnum, recdate, comment, vtype,
            naxis, dims, dbeg, dend):
        self.recnum = recnum
        self.recdate = recdate
        self.comment = comment
        self.vtype = vtype
        self.naxis = naxis
        self.dims = dims
        self.dbeg = dbeg
        self.dend = dend

    def __repr__(self):
        """Print contents of record header to terminal.
        """
        print('recnum={}, recdate={}, vtype={}, naxis={}, dims={}'.
                format(self.recnum, self.recdate, self.vtype,
                self.naxis, self.dims))
        print('dbeg={}, dend={}, comment="{}"'.
                format(self.dbeg, self.dend, self.comment))

class Dskfile:
    """Contents of a file in 'dsk' format.
    """

    def __init__(self, filename):
        self.filename = filename
        with open(self.filename, mode='rb') as file:
            self.bytes = file.read()
        self.bo = self.byteorder()
        self.reclist = self.reclist()

    def byteorder(self):
        """Use check bytes in the header to determine whether to swap bytes.
        The check bytes should yield 256 when read as an unsigned integer.
        """
        fmt = '>IH'
        hlen, check = struct.unpack(fmt, self.bytes[0:6])
        if check == 1:
            return '<'
        elif check == 256:
            return '>'
        else:
            raise Exception("{} not a 'dsk' file".format(self.filename))

    def rechead(self, recnum):
        """Read record header starting at the specified byte.
        Build a string with the date the record was written.
        Read the comment string from the header, if present.
        """
        ibeg, iend = self.reclist[recnum]
        fmt = self.bo + 'HBBBBBBH'
        (check, year, month, day, hours, mins, secs, clen) = \
                struct.unpack(fmt, self.bytes[ibeg:ibeg+10])
        ibeg += 10
        if year > 90:
            year += 1900
        else:
            year += 2000
        recdate = '{:04d}-{:02d}-{:02d}T{:02d}:{:02d}:{:02d}'. \
                format(year, month, day, hours, mins, secs)
        if clen == 0:
            comment = ''
        else:
            fmt = self.bo + '{}s'.format(clen)
            comment, = struct.unpack(fmt, self.bytes[ibeg:ibeg+clen])
            ibeg += clen
        fmt = self.bo + 'BB'
        (vtype, naxis) = struct.unpack(fmt, self.bytes[ibeg:ibeg+2])
        ibeg += 2
        fmt = self.bo + 'I' * naxis
        dims = struct.unpack(fmt, self.bytes[ibeg:ibeg+4*naxis])
        ibeg += 4 * naxis
        return Rechead(recnum, recdate, comment, vtype, naxis, dims, ibeg, iend)

    def recdata(self, rechead):
        """Read record data corresponding to the specified record header.
        Strings (vtype 7) are returned as a list of byte-array strings.
        Numeric data (vtypes 2, 3, 4, 5) are returned as an numpy array.
        Byteswapping inferred from the record header is applied to the data.
        Array dimensions are reversed from IDL standard to python standard.
        """
        dbeg = rechead.dbeg
        dend = rechead.dend
        vtype = rechead.vtype
        if vtype == 7:
            fmt = self.bo + 'h'
            recdata = list()
            i = dbeg
            while i+2 < dend:
                len, = struct.unpack(fmt, self.bytes[i:i+2])
                recdata.append(self.bytes[i+2:i+len+2])
                i += len + 2
            return recdata
        elif vtype == 2:
            dtype = 'u2'
        elif vtype == 3:
            dtype = 'u4'
        elif vtype == 4:
            dtype = 'f4'
        elif vtype == 5:
            dtype = 'f8'
        else:
            raise ValueError('unknown IDL vtype: {}'.format(vtype))
        bodt = self.bo + dtype
        dlen = (dend - dbeg) // np.dtype(bodt).itemsize
        return np.frombuffer(self.bytes, bodt, dlen, dbeg). \
                reshape(rechead.dims[::-1])

    def reclist(self):
        """Determine the first and last byte of each record, excluding the
        fortran-style record size at the beginning and end of each record.
        """
        fmt = self.bo + 'I'
        reclist = list()
        recnum = 0
        ibeg = 0
        while True:
            try:
                hlen, = struct.unpack(fmt, self.bytes[ibeg:ibeg+4])
            except struct.error:
                return reclist
            iend = ibeg + hlen + 4
            tlen, = struct.unpack(fmt, self.bytes[iend:iend+4])
            if hlen != tlen:
                raise Exception('record {}: hlen ({}) != tlen ({})'.
                        format(recnum, hlen, tlen))
            reclist.append((ibeg+4, iend))
            recnum += 1
            ibeg = iend + 4

    def record(self, recnum):
        """Extract contents of the specified record.
        """
        rechead = self.rechead(recnum)
        recdata = self.recdata(rechead)
        return recdata

def main():
    """Read the specified file in 'dsk' format and write 'fits' equivalent.
    Assume two records containing a reduced echelle spectrum.
    Assume record 0 contains primary data, returned as a numpy array.
    Assume record 1 contains FITS header, returned as a list of strings.
    """
    parser = argparse.ArgumentParser(
        description="Read data from a data file in 'dsk' format.",
        epilog='example: rdsk.py raab.51')
    parser.add_argument('filename', help="name of file in 'dsk' format")
    args = parser.parse_args()
    d = Dskfile(args.filename)
    fitsdata = d.record(0)
    fitshead = makehead(d.record(1))
    hdu = fits.PrimaryHDU(fitsdata, fitshead)
    hdulist = fits.HDUList([hdu])
    hdulist.writeto(args.filename+'.fits', clobber=True)

if __name__ == "__main__":
    main()
