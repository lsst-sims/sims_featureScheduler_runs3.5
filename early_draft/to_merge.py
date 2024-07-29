import numpy as np
import healpy as hp
import warnings

from rubin_scheduler.scheduler.utils import (
    EuclidOverlapFootprint,
    Footprint,
    Footprints,
    StepSlopes,
    set_default_nside,
)
from rubin_scheduler.utils import _hpid2_ra_dec, angular_separation
from astropy import units as u
from astropy.coordinates import SkyCoord


# Need to update LMC and SMC radius
class CustomAreaMap(EuclidOverlapFootprint):

    def __init__(
        self,
        nside=32,
        dust_limit=0.199,
        smoothing_cutoff=0.45,
        smoothing_beam=10,
        lmc_ra=80.893860,
        lmc_dec=-69.756126,
        lmc_radius=6,
        smc_ra=13.186588,
        smc_dec=-72.828599,
        smc_radius=4,
        scp_dec_max=-60,
        gal_long1=335,
        gal_long2=25,
        gal_lat_width_max=23,
        center_width=12,
        end_width=4,
        gal_dec_max=12,
        low_dust_dec_min=-70,
        low_dust_dec_max=15,
        adjust_halves=12,
        dusty_dec_min=-90,
        dusty_dec_max=15,
        eclat_min=-10,
        eclat_max=10,
        eclip_dec_min=0,
        nes_glon_limit=45.0,
        virgo_ra=186.75,
        virgo_dec=12.717,
        virgo_radius=8.75,
        euclid_contour_file=None,
    ):
        self.nside = nside
        self.hpid = np.arange(0, hp.nside2npix(nside))
        self.read_dustmap()

        self.lmc_ra = lmc_ra
        self.lmc_dec = lmc_dec
        self.lmc_radius = lmc_radius
        self.smc_ra = smc_ra
        self.smc_dec = smc_dec
        self.smc_radius = smc_radius

        self.virgo_ra = virgo_ra
        self.virgo_dec = virgo_dec
        self.virgo_radius = virgo_radius

        self.scp_dec_max = scp_dec_max

        self.gal_long1 = gal_long1
        self.gal_long2 = gal_long2
        self.gal_lat_width_max = gal_lat_width_max
        self.center_width = center_width
        self.end_width = end_width
        self.gal_dec_max = gal_dec_max

        self.low_dust_dec_min = low_dust_dec_min
        self.low_dust_dec_max = low_dust_dec_max
        self.adjust_halves = adjust_halves

        self.dusty_dec_min = dusty_dec_min
        self.dusty_dec_max = dusty_dec_max

        self.eclat_min = eclat_min
        self.eclat_max = eclat_max
        self.eclip_dec_min = eclip_dec_min
        self.nes_glon_limit = nes_glon_limit

        # Ra/dec in degrees and other coordinates
        self.ra, self.dec = hp.pix2ang(nside, self.hpid, lonlat=True)
        self.coord = SkyCoord(ra=self.ra * u.deg, dec=self.dec * u.deg, frame="icrs")
        self.eclip_lat = self.coord.barycentrictrueecliptic.lat.deg
        self.eclip_lon = self.coord.barycentrictrueecliptic.lon.deg
        self.gal_lon = self.coord.galactic.l.deg
        self.gal_lat = self.coord.galactic.b.deg

        # Set the low extinction area
        self.low_dust = np.where((self.dustmap < dust_limit), 1, 0)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", category=UserWarning)
            self.low_dust = hp.smoothing(self.low_dust, fwhm=np.radians(smoothing_beam))
        self.low_dust = np.where(self.low_dust > smoothing_cutoff, 1, 0)

        self.euclid_contour_file = euclid_contour_file

    def add_bulgy(self, filter_ratios, label="bulgy"):
        """Define a bulge region, where the 'bulge' is a series of
        circles set by points defined to match as best as possible the
        map requested by the SMWLV working group on galactic plane coverage.
        Implemented in v3.0.
        Updates self.healmaps and self.pix_labels.

        Parameters
        ----------
        filter_ratios : `dict` {`str`: `float`}
            Dictionary of weights per filter for the footprint.
        label : `str`, optional
            Label to apply to the resulting footprint
        """
        # Some RA, dec, radius points that
        # seem to cover the areas that are desired
        points = [
            [100.90, 9.55, 3],
            [84.92, -5.71, 3],
            [266.3, -29, 17],
            [279, -13, 10],
            [256, -45, 11],
            [155, -56.5, 6.5],
            [172, -62, 5],
            [190, -65, 5],
            [210, -64, 5],
            [242, -58, 6.5],
            [225, -60, 6.5],
        ]
        for point in points:
            dist = angular_separation(self.ra, self.dec, point[0], point[1])
            # Only change pixels where the label isn't already set.
            indx = np.where((dist < point[2]) & (self.pix_labels == ""))
            self.pix_labels[indx] = label
            for filtername in filter_ratios:
                self.healmaps[filtername][indx] = filter_ratios[filtername]

    def return_maps(
        self,
        magellenic_clouds_ratios={
            "u": 0.75,
            "g": 1.0,
            "r": 1.1,
            "i": 1.0,
            "z": 0.34,
            "y": 0.35,
        },
        scp_ratios={"u": 0.1, "g": 0.175, "r": 0.1, "i": 0.135, "z": 0.046, "y": 0.047},
        nes_ratios={"g": 0.255, "r": 0.33, "i": 0.33, "z": 0.23},
        dusty_plane_ratios={
            "u": 0.093,
            "g": 0.26,
            "r": 0.26,
            "i": 0.26,
            "z": 0.26,
            "y": 0.093,
        },
        low_dust_ratios={"u": 0.35, "g": 0.4, "r": 1.0, "i": 1.0, "z": 0.9, "y": 0.9},
        bulge_ratios={"u": 0.177, "g": 0.98, "r": 1.03, "i": 1.03, "z": 0.98, "y": 0.226},
        virgo_ratios={"u": 0.35, "g": 0.4, "r": 1.0, "i": 1.0, "z": 0.9, "y": 0.9},
        euclid_ratios={"u": 0.35, "g": 0.4, "r": 1.0, "i": 1.0, "z": 0.9, "y": 0.9},
    ):

        # Array to hold the labels for each pixel
        self.pix_labels = np.zeros(hp.nside2npix(self.nside), dtype="U20")
        self.healmaps = np.zeros(
            hp.nside2npix(self.nside),
            dtype=list(zip(["u", "g", "r", "i", "z", "y"], [float] * 7)),
        )

        # Note, order here matters.
        # Once a HEALpix is set and labled, subsequent add_ methods
        # will not override that pixel.
        self.add_magellanic_clouds(magellenic_clouds_ratios)
        self.add_lowdust_wfd(low_dust_ratios)
        self.add_virgo_cluster(virgo_ratios)
        self.add_bulgy(bulge_ratios)
        self.add_nes(nes_ratios)
        self.add_dusty_plane(dusty_plane_ratios)
        self.add_euclid_overlap(euclid_ratios)
        self.add_scp(scp_ratios)

        return self.healmaps, self.pix_labels


def make_rolling_footprints(
    fp_hp=None,
    mjd_start=60218.0,
    sun_ra_start=3.27717639,
    nslice=2,
    scale=0.8,
    nside=32,
    wfd_indx=None,
    order_roll=0,
    n_cycles=None,
    n_constant_start=2,
    n_constant_end=6,
    verbose=False,
    uniform=True,
):
    """
    Generate rolling footprints

    Parameters
    ----------
    fp_hp : dict-like
        A dict with filtername keys and HEALpix map values
    mjd_start : `float`
        The starting date of the survey.
    sun_ra_start : `float`
        The RA of the sun at the start of the survey
    nslice : `int`
        How much to slice the sky up. Can be 2, 3, 4, or 6.
    scale : `float`
        The strength of the rolling, value of 1 is full power rolling.
        Zero is no rolling.
    wfd_indx : array of ints
        The indices of the HEALpix map that are to be included in the rolling.
    order_roll : `int`
        Change the order of when bands roll. Default 0.
    n_cycles : `int`
        Number of complete rolling cycles to attempt. If None, defaults to 3
        full cycles for nslice=2, 2 cycles for nslice=3 or 4, and 1 cycle for
        nslice=6.
    n_constant_start : `int`
        The number of constant non-rolling seasons to start with. Anything
        less than 2 will start rolling too early near Y1. Defaults to 2.
    n_constant_end : `int`
        The number of constant seasons to end the survey with. Defaults to 6.

    Returns
    -------
    Footprints object
    """

    nc_default = {2: 3, 3: 2, 4: 2, 6: 1}
    if n_cycles is None:
        n_cycles = nc_default[nslice]

    hp_footprints = fp_hp

    D = 1.0 - scale
    U = nslice - D * (nslice - 1)

    start = [1.0] * n_constant_start
    # After n_cycles, just go to no-rolling for 6 years.
    end = [1.0] * n_constant_end

    rolling = [U] + [D] * (nslice - 1)
    rolling = np.roll(rolling, order_roll).tolist()

    all_slopes = []
    if uniform:
        extra_constant = [1]
    else:
        extra_constant = []

    for i in range(nslice):
        _roll = np.roll(rolling, i).tolist() + extra_constant
        all_slopes.append(
            start + _roll * n_cycles + end
        )
    for i in range(nslice):
        _roll = np.roll(rolling, i).tolist() + extra_constant
        _roll = [_roll[-1]] + _roll[1:-1] + [_roll[0]]
        all_slopes.append(
            start + _roll * n_cycles + end
            )
    dvals = {
        1: "1",
        D: "D",
        U: "U",
    }

    abc = ["a", "b", "c", "d", "e", "f", "g", "h"]
    slice_names = ["slice %s" % abc[i] for i in range(nslice)]
    for i, s in enumerate(all_slopes):
        if i >= nslice:
            sname = (
                slice_names[i-nslice]
                + " w/ ra - sun_ra in [90, 270]"
            )
        else:
            sname = (
                slice_names[i]
                + " w/ ra - sun_ra in [270, 90]"
            )
        if verbose:
            print(sname + ": " + " ".join([dvals[x] for x in s]))

    fp_non_wfd = Footprint(mjd_start, sun_ra_start=sun_ra_start, nside=nside)
    rolling_footprints = []
    for i in range(len(all_slopes)):
        step_func = StepSlopes(rise=all_slopes[i])
        rolling_footprints.append(
            Footprint(
                mjd_start,
                sun_ra_start=sun_ra_start,
                step_func=step_func,
                nside=nside,
            )
        )

    wfd = hp_footprints["r"] * 0
    if wfd_indx is None:
        wfd_indx = np.where(hp_footprints["r"] == 1)[0]

    wfd[wfd_indx] = 1
    non_wfd_indx = np.where(wfd == 0)[0]

    if uniform:
        split_wfd_indices = slice_quad_galactic_cut(
            hp_footprints, nslice=nslice, wfd_indx=wfd_indx,
            ra_range=(sun_ra_start + 1.5 * np.pi, sun_ra_start + np.pi/2),
        )

        split_wfd_indices_delayed = slice_quad_galactic_cut(
            hp_footprints, nslice=nslice, wfd_indx=wfd_indx,
            ra_range=(sun_ra_start + np.pi / 2, sun_ra_start + 1.5 * np.pi),
        )
    else:
        split_wfd_indices = slice_quad_galactic_cut(hp_footprints, nslice=nslice, wfd_indx=wfd_indx)

    for key in hp_footprints:
        temp = hp_footprints[key] + 0
        temp[wfd_indx] = 0
        fp_non_wfd.set_footprint(key, temp)

        for i in range(nslice):
            # make a copy of the current filter
            temp = hp_footprints[key] + 0
            # Set the non-rolling area to zero
            temp[non_wfd_indx] = 0

            indx = split_wfd_indices[i]
            # invert the indices
            ze = temp * 0
            ze[indx] = 1
            temp = temp * ze
            rolling_footprints[i].set_footprint(key, temp)
        if uniform:
            for _i in range(nslice, nslice*2):
                # make a copy of the current filter
                temp = hp_footprints[key] + 0
                # Set the non-rolling area to zero
                temp[non_wfd_indx] = 0

                indx = split_wfd_indices_delayed[_i-nslice]
                # invert the indices
                ze = temp * 0
                ze[indx] = 1
                temp = temp * ze
                rolling_footprints[_i].set_footprint(key, temp)

    result = Footprints([fp_non_wfd] + rolling_footprints)
    return result


def _is_in_ra_range(ra, low, high):
    _low = low % (2.0 * np.pi)
    _high = high % (2.0 * np.pi)
    if _low <= _high:
        return (ra >= _low) & (ra <= _high)
    else:
        return (ra >= _low) | (ra <= _high)


def slice_quad_galactic_cut(
    target_map, nslice=2, wfd_indx=None, ra_range=None
):
    """
    Helper function for generating rolling footprints

    Parameters
    ----------
    target_map : dict of HEALpix maps
        The final desired footprint as HEALpix maps. Keys are filter names
    nslice : `int`
        The number of slices to make, can be 2 or 3.
    wfd_indx : array of ints
        The indices of target_map that should be used for rolling.
        If None, assumes the rolling area should be where target_map['r'] == 1.
    ra_range : tuple of floats, optional
        If not None, then the indices are restricted to the given RA range
        in radians.
    """

    ra, dec = ra_dec_hp_map(nside=hp.npix2nside(target_map["r"].size))

    coord = SkyCoord(ra=ra * u.rad, dec=dec * u.rad)
    _, gal_lat = coord.galactic.l.deg, coord.galactic.b.deg

    indx_north = np.intersect1d(np.where(gal_lat >= 0)[0], wfd_indx)
    indx_south = np.intersect1d(np.where(gal_lat < 0)[0], wfd_indx)

    splits_north = slice_wfd_area_quad(
        target_map, nslice=nslice, wfd_indx=indx_north)
    splits_south = slice_wfd_area_quad(
        target_map, nslice=nslice, wfd_indx=indx_south)

    slice_indx = []
    for j in np.arange(nslice):
        indx_temp = []
        for i in np.arange(j + 1, nslice * 2 + 1, nslice):
            indx_temp += indx_north[
                splits_north[i - 1]: splits_north[i]
            ].tolist()
            indx_temp += indx_south[
                splits_south[i - 1]: splits_south[i]
            ].tolist()
        slice_indx.append(indx_temp)

    if ra_range is not None:
        ra_indx = np.where(_is_in_ra_range(ra, *ra_range))[0]
        for j in range(nslice):
            slice_indx[j] = np.intersect1d(ra_indx, slice_indx[j])

    return slice_indx


def slice_wfd_area_quad(target_map, nslice=2, wfd_indx=None):
    """
    Divide a healpix map in an intelligent way

    Parameters
    ----------
    target_map : dict of HEALpix arrays
        The input map to slice
    nslice : int
        The number of slices to divide the sky into (gets doubled).
    wfd_indx : array of int
        The indices of the healpix map to consider as part of the WFD area
        that will be split.
        If set to None, the pixels where target_map['r'] == 1 are
        considered as WFD.
    """
    nslice2 = nslice * 2

    wfd = target_map["r"] * 0
    if wfd_indx is None:
        wfd_indices = np.where(target_map["r"] == 1)[0]
    else:
        wfd_indices = wfd_indx
    wfd[wfd_indices] = 1
    wfd_accum = np.cumsum(wfd)
    split_wfd_indices = np.floor(
        np.max(wfd_accum) / nslice2 * (np.arange(nslice2) + 1)
    ).astype(int)
    split_wfd_indices = split_wfd_indices.tolist()
    split_wfd_indices = [0] + split_wfd_indices

    return split_wfd_indices


def ra_dec_hp_map(nside=None):
    """
    Return all the RA,dec points for the centers of a healpix map, in radians.
    """
    if nside is None:
        nside = set_default_nside()
    ra, dec = _hpid2_ra_dec(nside, np.arange(hp.nside2npix(nside)))
    return ra, dec
