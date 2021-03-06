"""
Contains triallelic Spectrum object
"""
import os
import numpy as np
import dadi
import numerics

class TLSpectrum(np.ma.masked_array):
    """
    Represents a two-locus frequency spectrum.

    The constructor has the format:
        fs = dadi.Triallele.TLSpectrum(data, mask, mask_infeasible,
                                        data_folded,
                                        extrap_x, extrap_t)

        data: The triallelic frequency spectrum data
        mask: An optional array of the same size as data, similar to dadi.Spectrum
        data_folded: If True, it is assumed that the input data is folded
        check_folding: If True and data_folded=True, the data and
                       mask will be checked to ensure they are consistent
        extrap_x: Optional floating point value specifying x value to use in
                  extrapolation.
        extrap_t: Optional floating point value specifying t value to use in
                  extrapolation.
    """
    def __new__(subtype, data, mask=np.ma.nomask, mask_infeasible=True, 
                data_folded=None, check_folding=True,
                dtype=float, copy=True, fill_value=np.nan, keep_mask=True,
                shrink=True, extrap_x=None, extrap_t=None):
        data = np.asanyarray(data)
        
        if mask is np.ma.nomask:
            mask = np.ma.make_mask_none(data.shape)
        
        subarr = np.ma.masked_array(data, mask=mask, dtype=dtype, copy=copy,
                                       fill_value=fill_value, keep_mask=True, 
                                       shrink=True)
        subarr = subarr.view(subtype)
        if hasattr(data, 'folded'):
            if data_folded is None or data_folded == data.folded:
                subarr.folded = data.folded
            elif data_folded != data.folded:
                raise ValueError('Data does not have same folding status as '
                                 'was called for in TLSpectrum constructor.')
        elif data_folded is not None:
            subarr.folded = data_folded
        else:
            subarr.folded = False
        
        ### XXX To do: ensure that all goes well when creating the TLSpectrum, come
        ###     back to this
        # Check that if we're declaring that the input data is folded, it actually is,
        # and the mask reflects this.        
                
        if mask_infeasible:
            subarr.mask_infeasible()
        
        subarr.extrap_x = extrap_x
        subarr.extrap_t = extrap_t
        
        return subarr
        
    # See https://scipy.github.io/old-wiki/pages/Subclasses.html for information on
    # __array_finalize__ and __array_wrap__ methods.
    #
    # We need these methods to ensure extra attributes get copied along when
    # we do arithmetic on the FS.
    def __array_finalize__(self, obj):
        if obj is None: 
            return
        np.ma.masked_array.__array_finalize__(self, obj)
        self.folded = getattr(obj, 'folded', 'unspecified')
        self.extrap_x = getattr(obj, 'extrap_x', None)
        self.extrap_t = getattr(obj, 'extrap_t', None)
    def __array_wrap__(self, obj, context=None):
        result = obj.view(type(self))
        result = np.ma.masked_array.__array_wrap__(self, obj, 
                                                      context=context)
        result.folded = self.folded
        result.extrap_t = self.extrap_t
        return result
    def _update_from(self, obj):
        np.ma.masked_array._update_from(self, obj)
        if hasattr(obj, 'folded'):
            self.folded = obj.folded
        if hasattr(obj, 'extrap_x'):
            self.extrap_x = obj.extrap_x
        if hasattr(obj, 'extrap_t'):
            self.extrap_t = obj.extrap_t
    # masked_array has priority 15.
    __array_priority__ = 20

    def __repr__(self):
        return 'TLSpectrum(%s, folded=%s)'\
                % (str(self), str(self.folded))

    def mask_infeasible(self):
        """
        Mask any infeasible entries.
        """
        ns = len(self)-1
        self.mask[0,0,0] = True
        self.mask[0,:,0] = True
        self.mask[0,0,:] = True
        for ii in range(len(self)):
            for jj in range(len(self)):
                for kk in range(len(self)):
                    if ii+jj+kk > ns:
                        self.mask[ii,jj,kk] = True

        for ii in range(len(self)):
            self.mask[ii,ns-ii,0] = True
            self.mask[ii,0,ns-ii] = True
        
        return self
    
    def unfold(self):
        if not self.folded:
            raise ValueError('Input Spectrum is not folded.')
        data = self.data
        unfolded = TLSpectrum(data, mask_infeasible=True)
        unfolded.extrap_x = self.extrap_x
        unfolded.extrap_t = self.extrap_t
        return unfolded

    def _get_sample_size(self):
        return np.asarray(self.shape)[0] - 1
    sample_size = property(_get_sample_size)
    
    def _ensure_shape_and_dimension(self):
        """
        Ensure that fs has Npop dimensions.
        """
        pass
    
    # Make from_file a static method, so we can use it without an instance.
    @staticmethod
    def from_file(fid, mask_infeasible=True, return_comments=False):
        """
        Read frequency spectrum from file.

        fid: string with file name to read from or an open file object.
        mask_infeasible: If True, mask the infeasible entries in the triallelic spectrum.
        return_comments: If true, the return value is (fs, comments), where
                         comments is a list of strings containing the comments
                         from the file (without #'s).

        See to_file method for details on the file format.
        """
        newfile = False
        # Try to read from fid. If we can't, assume it's something that we can
        # use to open a file.
        if not hasattr(fid, 'read'):
            newfile = True
            fid = file(fid, 'r')

        line = fid.readline()
        # Strip out the comments
        comments = []
        while line.startswith('#'):
            comments.append(line[1:].strip())
            line = fid.readline()

        # Read the shape of the data
        shape,folded,extrap_x,extrap_t = line.split()
        shape = [int(shape)+1,int(shape)+1,int(shape)+1]

        data = np.fromstring(fid.readline().strip(), 
                                count=np.product(shape), sep=' ')
        # fromfile returns a 1-d array. Reshape it to the proper form.
        data = data.reshape(*shape)

        maskline = fid.readline().strip()
        mask = np.fromstring(maskline, 
                                count=np.product(shape), sep=' ')
        mask = mask.reshape(*shape)
        
        if folded == 'folded':
            folded = True
        else:
            folded = False
        if extrap_x == 'None':
            extrap_x = None
        else:
            extrap_x = float(extrap_x)
        if extrap_t == 'None':
            extrap_t = None
        else:
            extrap_t = float(extrap_t)

        # If we opened a new file, clean it up.
        if newfile:
            fid.close()

        fs = TLSpectrum(data, mask, mask_infeasible, data_folded=folded)
        fs.extrap_x = extrap_x
        fs.extrap_t = extrap_t
        if not return_comments:
            return fs
        else:
            return fs,comments
    
    def to_file(self, fid, precision=16, comment_lines=[], foldmaskinfo=True, extrapinfo=True):
        """
        Write frequency spectrum to file.
    
        fid: string with file name to write to or an open file object.
        precision: precision with which to write out entries of the SFS. (They 
                   are formated via %.<p>g, where <p> is the precision.)
        comment lines: list of strings to be used as comment lines in the header
                       of the output file.
        foldmaskinfo: If False, folding and mask and population label
                      information will not be saved. This conforms to the file
                      format for dadi versions prior to 1.3.0.

        The file format is:
            # Any number of comment lines beginning with a '#'
            A single line containing N integers giving the dimensions of the fs
              array. So this line would be '5 5 3' for an SFS that was 5x5x3.
              (That would be 4x4x2 *samples*.)
            On the *same line*, the string 'folded' or 'unfolded' 
              denoting the folding status of the array
            A single line giving the array elements. The order of elements is 
              e.g.: fs[0,0,0] fs[0,0,1] fs[0,0,2] ... fs[0,1,0] fs[0,1,1] ...
            A single line giving the elements of the mask in the same order as
              the data line. '1' indicates masked, '0' indicates unmasked.
        """
        # Open the file object.
        newfile = False
        if not hasattr(fid, 'write'):
            newfile = True
            fid = file(fid, 'w')

        # Write comments
        for line in comment_lines:
            fid.write('# ')
            fid.write(line.strip())
            fid.write(os.linesep)

        # Write out the shape of the fs
        fid.write('{0} '.format(self.sample_size))

        if foldmaskinfo:
            if not self.folded:
                fid.write('unfolded ')
            else:
                fid.write('folded ')
        
        if extrapinfo:
            if not self.extrap_x:
                fid.write('None ')
            else:
                fid.write('{0} '.format(self.extrap_x))
            if not self.extrap_t:
                fid.write('None')
            else:
                fid.write('{0}'.format(self.extrap_t))

        fid.write(os.linesep)

        # Write the data to the file
        self.data.tofile(fid, ' ', '%%.%ig' % precision)
        fid.write(os.linesep)

        if foldmaskinfo:
            # Write the mask to the file
            np.asarray(self.mask,int).tofile(fid, ' ')
            fid.write(os.linesep)

        # Close file
        if newfile:
            fid.close()

    tofile = to_file

    def marginalA(self):
        """
        Marginal 1D frequency spectrum for A locus.
        """
        ns = self.shape[0] - 1
        marg = dadi.Spectrum(np.zeros(ns+1))
        for fAB in range(ns):
            for fAb in range(ns-fAB):
                marg[fAB+fAb] += self[fAB,fAb,:].sum()

        marg.extrap_x = self.extrap_x
        marg.extrap_t = self.extrap_t
        return marg

    def marginalB(self):
        """
        Marginal 1D frequency spectrum for B locus.
        """
        ns = self.shape[0] - 1
        marg = dadi.Spectrum(np.zeros(ns+1))
        for fAB in range(ns):
            for faB in range(ns-fAB):
                marg[fAB+faB] += self[fAB,:,faB].sum()

        marg.extrap_x = self.extrap_x
        marg.extrap_t = self.extrap_t
        return marg

    def mean_r2(self):
        """
        Mean of normalized squared correlation coefficient between A and B loci.
        """
        ns = self.shape[0] - 1
        norm = self.sum()
        Dbin, r2bin = numerics.LD_per_bin(ns)
        return (self*r2bin).sum()/self.sum()

    def fold(self):
        if self.folded:
            raise ValueError('Input Spectrum is already folded.')
        ns = self.shape[0] - 1
        folded = 0*self
        for ii in range(ns+1):
            for jj in range(ns+1):
                for kk in range(ns+1):
                    if self.mask[ii,jj,kk]:
                        continue
                    p = ii + jj
                    q = ii + kk
                    if p > ns/2 and q > ns/2:
                        # Switch A/a and B/b, so AB becomes ab, Ab becomes aB, etc
                        folded[ns-ii-jj-kk,kk,jj] = self[ns-ii-jj-kk,kk,jj] + self[ii,jj,kk]
                        folded.mask[ii,jj,kk] = True
                    elif p > ns/2:
                        # Switch A/a, so AB -> aB, Ab -> ab, aB -> AB, and ab -> Ab
                        folded[kk,ns-ii-jj-kk,ii] = self[kk,ns-ii-jj-kk,ii] + self[ii,jj,kk]
                        folded.mask[ii,jj,kk] = True
                    elif q > ns/2:
                        # Switch B/b, so AB -> Ab, Ab -> AB, aB -> ab, and ab -> aB
                        folded[jj,ii,ns-ii-jj-kk] = self[jj,ii,ns-ii-jj-kk] + self[ii,jj,kk]
                        folded.mask[ii,jj,kk] = True
        
        folded.folded = True
        folded.extrap_x = self.extrap_x
        folded.extrap_t = self.extrap_t
        return folded
    
    # Ensures that when arithmetic is done with TLSpectrum objects,
    # attributes are preserved. For details, see similar code in
    # dadi.Spectrum_mod
    for method in ['__add__','__radd__','__sub__','__rsub__','__mul__',
                   '__rmul__','__div__','__rdiv__','__truediv__','__rtruediv__',
                   '__floordiv__','__rfloordiv__','__rpow__','__pow__']:
        exec("""
def %(method)s(self, other):
    self._check_other_folding(other)
    if isinstance(other, np.ma.masked_array):
        newdata = self.data.%(method)s (other.data)
        newmask = np.ma.mask_or(self.mask, other.mask)
    else:
        newdata = self.data.%(method)s (other)
        newmask = self.mask
    if hasattr(other, 'extrap_x') and self.extrap_x != other.extrap_x:
        extrap_x = None
    else:
        extrap_x = self.extrap_x
    if hasattr(other, 'extrap_t') and self.extrap_t != other.extrap_t:
        extrap_t = None
    else:
        extrap_t = self.extrap_t
    outfs = self.__class__.__new__(self.__class__, newdata, newmask, 
                                   mask_infeasible=False, 
                                   data_folded=self.folded,
                                   extrap_x=extrap_x, extrap_t=extrap_t)
    return outfs
""" % {'method':method})

    # Methods that modify the Spectrum in-place.
    for method in ['__iadd__','__isub__','__imul__','__idiv__',
                   '__itruediv__','__ifloordiv__','__ipow__']:
        exec("""
def %(method)s(self, other):
    self._check_other_folding(other)
    if isinstance(other, np.ma.masked_array):
        self.data.%(method)s (other.data)
        self.mask = np.ma.mask_or(self.mask, other.mask)
    else:
        self.data.%(method)s (other)
    if hasattr(other, 'extrap_x') and self.extrap_x != other.extrap_x:
        self.extrap_x = None
    if hasattr(other, 'extrap_t') and self.extrap_t != other.extrap_t:
        self.extrap_t = None
    return self
""" % {'method':method})

    def _check_other_folding(self, other):
        """
        Ensure other Spectrum has same .folded status
        """
        if isinstance(other, self.__class__)\
           and (other.folded != self.folded):
            raise ValueError('Cannot operate with a folded Spectrum and an '
                             'unfolded one.')


# Allow TLSpectrum objects to be pickled. 
# See http://effbot.org/librarybook/copy-reg.htm
import copy_reg
def TLSpectrum_pickler(fs):
    # Collect all the info necessary to save the state of a TLSpectrum
    return TLSpectrum_unpickler, (fs.data, fs.mask, fs.folded,
                                   fs.extrap_x, fs.extrap_t)
def TLSpectrum_unpickler(data, mask, folded,
                          extrap_x, extrap_t):
    # Use that info to recreate the TLSpectrum
    return TLSpectrum(data, mask, mask_infeasible=False,
                       data_folded=folded,
                       extrap_x=extrap_x, extrap_t=extrap_t)
copy_reg.pickle(TLSpectrum, TLSpectrum_pickler, TLSpectrum_unpickler)
