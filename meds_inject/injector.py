import numpy as np
import galsim
import ngmix
import fitsio
from .pbar import prange

SCALE = 0.263
PSF_NOISE = 1.0e-6


def make_obj(args):
    """
    make the galsim object to be injected
    """
    flux = 10.0**(0.4*(30.0 - args.mag))
    return galsim.Gaussian(
        fwhm=args.fwhm,
        flux=flux,
    ).shear(
        g1=args.g1,
        g2=args.g2,
    )


def make_psf(args):
    """
    make the psf to be used for convolutions
    """
    return galsim.Gaussian(
        fwhm=args.psf_fwhm,
    )


def make_image(obj, dims, offset=None):
    """
    draw an image, possibly with an offset
    """
    return obj.drawImage(
        nx=dims[1],
        ny=dims[0],
        scale=SCALE,
        offset=offset,
    ).array


def add_noise(*, rng, image, noise):
    """
    add gaussian noise to an image
    """
    image += rng.normal(
        scale=noise,
        size=image.shape,
    )


def make_obj_image(*, args, rng, cat, iobj, icut):
    """
    make the image and reset the catalog values
    cutout_row, cutout_col
    """
    psf_obj = make_psf(args)
    obj0 = make_obj(args)
    obj = galsim.Convolve(obj0, psf_obj)

    offset = rng.uniform(low=-0.5, high=0.5, size=2)

    box_size = cat['box_size'][iobj]

    im = make_image(obj, [box_size]*2, offset=offset)

    cen = (np.array(im.shape)-1.0)/2.0 + offset

    add_noise(rng=rng, image=im, noise=args.noise)

    cat['cutout_row'][iobj, icut] = cen[0]
    cat['cutout_col'][iobj, icut] = cen[1]

    return im.astype('f4')


def make_psf_image(*, args, rng, cat, iobj, icut):
    """
    make the image and reset the catalog values
    psf_cutout_row, psf_cutout_col
    """
    psf_obj = make_psf(args)

    row_size = cat['psf_row_size'][iobj, icut]
    col_size = cat['psf_col_size'][iobj, icut]

    dims = [row_size, col_size]

    im = make_image(psf_obj, dims)

    cen = (np.array(im.shape)-1.0)/2.0

    add_noise(rng=rng, image=im, noise=PSF_NOISE)

    cat['psf_cutout_row'][iobj, icut] = cen[0]
    cat['psf_cutout_col'][iobj, icut] = cen[1]

    return im.astype('f4')


def make_weight_image(*, args, rng, cat, iobj):
    """
    make the weight image
    """

    box_size = cat['box_size'][iobj]
    weight = np.zeros([box_size]*2, dtype='f4')
    weight[:, :] = 1.0/args.noise**2
    return weight


def make_bmask_image(*, args, rng, cat, iobj):
    """
    make the bmask image
    """

    box_size = cat['box_size'][iobj]
    return np.zeros([box_size]*2, dtype='i4')


def make_seg_image(*, args, rng, cat, iobj):
    """
    make the seg image, all pixels owned by
    the indicated object
    """

    box_size = cat['box_size'][iobj]
    seg = np.zeros([box_size]*2, dtype='i4')
    seg[:, :] = cat['number'][iobj]
    return seg


def print_progress(num, i):
    if i == 0 or (i + 1) % 100 == 0:
        print('%d/%d  %g%%' % (i+1, num, (i+1)/num*100))


def get_cutout(*, args, rng, cat, iobj, icut, cutout_type):
    """
    get a simulated cutout, possibly resetting values
    in the catalog such as cutout_row, cutout_col
    """
    if cutout_type == 'image_cutouts':
        image = make_obj_image(
            args=args,
            rng=rng,
            cat=cat,
            iobj=iobj,
            icut=icut,
        )
    elif cutout_type == 'weight_cutouts':
        image = make_weight_image(
            args=args,
            rng=rng,
            cat=cat,
            iobj=iobj,
        )
    elif cutout_type == 'bmask_cutouts':
        image = make_bmask_image(
            args=args,
            rng=rng,
            cat=cat,
            iobj=iobj,
        )
    elif cutout_type == 'seg_cutouts':
        image = make_seg_image(
            args=args,
            rng=rng,
            cat=cat,
            iobj=iobj,
        )
    else:
        raise ValueError('bad cutout_type: %s' % cutout_type)
    return image


def write_cutouts(*, args, fits, cat, cutout_type, rng):
    """
    generic writer for cutouts other than psf
    """
    print('writing', cutout_type)

    hdu = fits[cutout_type]

    if args.ntest is not None:
        num = args.ntest
    else:
        num = cat.size

    for iobj in prange(num):
        print_progress(cat.size, iobj)

        for icut in range(cat['ncutout'][iobj]):
            image = get_cutout(
                args=args,
                rng=rng,
                cat=cat,
                iobj=iobj,
                icut=icut,
                cutout_type=cutout_type,
            )
            start_row = cat['start_row'][iobj, icut]
            hdu.write(image, start=start_row)


def write_psf_cutouts(*, args, fits, cat, rng):
    """
    special writer for psfs, since they can be
    non-squar
    """
    print('writing psfs')

    hdu = fits['psf']

    if args.ntest is not None:
        num = args.ntest
    else:
        num = cat.size

    for iobj in prange(num):
        print_progress(cat.size, iobj)

        for icut in range(cat['ncutout'][iobj]):
            image = make_psf_image(
                args=args,
                rng=rng,
                cat=cat,
                iobj=iobj,
                icut=icut,
            )
            start_row = cat['psf_start_row'][iobj, icut]
            hdu.write(image, start=start_row)


def write_cat(*, fits, cat):
    """
    overwrite the catalog entries
    """
    hdu = fits['object_data']
    hdu.write(cat, start=0)


def set_jacobian(cat):
    """
    set all jacobian entries according to our
    chosen pixel scale
    """
    jac = ngmix.DiagonalJacobian(
        row=35, col=35, scale=SCALE,
    )

    cat['dudrow'] = jac.dudrow
    cat['dudcol'] = jac.dudcol
    cat['dvdrow'] = jac.dvdrow
    cat['dvdcol'] = jac.dvdcol


def inject(args):

    rng = np.random.RandomState(args.seed)

    with fitsio.FITS(args.file, 'rw') as fits:
        cat = fits['object_data'].read()

        set_jacobian(cat)

        cutout_types = [
            'image_cutouts',
            'weight_cutouts',
            'bmask_cutouts',
            'seg_cutouts',
            'psf',
        ]
        for cutout_type in cutout_types:
            if cutout_type == 'psf':
                write_psf_cutouts(
                    args=args,
                    fits=fits,
                    cat=cat,
                    rng=rng,
                )
            else:
                write_cutouts(
                    args=args,
                    fits=fits,
                    cat=cat,
                    cutout_type=cutout_type,
                    rng=rng,
                )

        write_cat(fits=fits, cat=cat)
