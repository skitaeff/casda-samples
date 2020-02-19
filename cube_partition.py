#!/usr/bin/env python3.7 -u

# The script orders CASDA to partition a cube without downloading
# Python 3.7
# Author: Slava Kitaeff
# Copywrite: AusSRC, 2019

from __future__ import print_function, division

import argparse
from configparser import ConfigParser
from typing import List

import casda
import configparser
from astropy.io.votable import parse
from astropy.io.votable import parse_single_table
import numpy as np
import time

def parseargs():
    """
    Parse the command line arguments
    :return: An args map with the parsed arguments
    """
    parser = argparse.ArgumentParser(description="Order CASDA cutouts without downloading")
    parser.add_argument("-u", "--opal_username", help="Your user name on the ATNF's online proposal system (normally an email address)")
    parser.add_argument("-p", "--opal_password", help="Your password on the ATNF's online proposal system")
    parser.add_argument("--password_file", help="The file holding your password for the ATNF's online proposal system")
    parser.add_argument("-c", "--conf_file", help="The configuration file holding parameters of partitioning")
    parser.add_argument("-f", "--cube_file_name", help="The image cube file name to be accessed in CASDA")
    parser.add_argument("-d", "--destination_dir", help="Local destination directory on processing node")
    parser.add_argument("--download", help="Download all cutouts and checksums")

    args = parser.parse_args()
    return args

def get_vo_table(file_name: str, username: str, password: str) -> object:

    query = "SELECT * from ivoa.obscore where filename = '" + file_name + "'"
    image_cube_votable = parse(casda.sync_tap_query_xml(query, username, password), pedantic=False)

    return image_cube_votable.get_table_by_id('results').array

def generate_cutouts_array(config_file: str, cube_vo_table: object) -> object:
    config: ConfigParser = configparser.ConfigParser()
    config.read_file(open(config_file))
    if 'Partitions' in config:
        RA_Partitions: int = int(config['Partitions']['RA_Partitions'])
        RA_Overlap: int = int(config['Partitions']['RA_Overlap'])
        Dec_Partitions: int = int(config['Partitions']['Dec_Partitions'])
        Dec_Overlap: int = int(config['Partitions']['Dec_Overlap'])
        Frequency_Partitions: int = int(config['Partitions']['Frequency_Partitions'])
        Frequency_Overlap: int = int(config['Partitions']['Frequency_Overlap'])
    else:
        raise RuntimeError('Could not find section Partitions in configuration file.')

    s_xel1: int = int(cube_vo_table['s_xel1']) #number of pixels on RA axes
    s_xel2: int = int(cube_vo_table['s_xel2']) #number of pixels on Dec axes
    em_xel: int = int(cube_vo_table['em_xel']) #number of frequency channels
    em_min: float = float(cube_vo_table['em_min']) #first channel in meters
    em_max: float = float(cube_vo_table['em_max']) #last channel in meters
    spectral_resolution: float = (em_max - em_min) / em_xel #spectral resolution per channel
    s_ra: float = float(cube_vo_table['s_ra']) #central RA in degrees
    s_dec: float = float(cube_vo_table['s_dec']) #central Dec in degrees
    spatial_resolution: float = float(cube_vo_table['s_resolution'])/360 #spacial resolution in degrees, s_resolution is in arcsec

# todo: it's unclear what units are actually used in spatial_resolution
    ra_1 = s_ra - spatial_resolution * s_xel1 / 2 # RA of the first pixel
    dec_1 = s_dec - spatial_resolution * s_xel2 / 2 # Dec of the first pixel

    dRA = s_xel1/RA_Partitions
    dDec = s_xel2/Dec_Partitions
    dFreq = em_xel/Frequency_Partitions

    ra_min_list = []
    ra_max_list = []
    # assign the first RA partition
    ra_min_list.append(1)
    ra_max_list.append(dRA + RA_Overlap)
    # assign the middle RA partition
    for i in range(1, RA_Partitions-1, 1):
        ra_min_list.append(dRA * i - RA_Overlap)
        ra_max_list.append(dRA * (i+1) + RA_Overlap)
    # assign the last RA partition
    ra_min_list.append(s_xel1 - (dRA + RA_Overlap))
    ra_max_list.append(s_xel1)
    # convert to degrees
    for i in range(0, RA_Partitions, 1):
        ra_min_list[i] = ra_1 + spatial_resolution * ra_min_list[i]
        ra_max_list[i] = ra_1 + spatial_resolution * ra_max_list[i]

    dec_min_list = []
    dec_max_list = []
    # assign the first Dec partition
    dec_min_list.append(1)
    dec_max_list.append(dDec + Dec_Overlap)
    # assign the middle Dec partition
    for i in range(1, Dec_Partitions-1, 1):
        dec_min_list.append(dDec * i - Dec_Overlap)
        dec_max_list.append(dDec * (i+1) + Dec_Overlap)
    # assign the last Dec partition
    dec_min_list.append(s_xel2 - (dDec + Dec_Overlap))
    dec_max_list.append(s_xel2)
    # convert to degrees
    for i in range(0, Dec_Partitions, 1):
        dec_min_list[i] = dec_1 + spatial_resolution * dec_min_list[i]
        dec_max_list[i] = dec_1 + spatial_resolution * dec_max_list[i]

    range_params = []
    for i in range(0, RA_Partitions, 1):
        for j in range(0, Dec_Partitions, 1):
            range_params.append("RANGE " + str(ra_min_list[i]) + " " + str(ra_max_list[i]) + " "+ str(dec_min_list[j]) + " "  + str(dec_max_list[j]))

    fre_min_list = []
    fre_max_list = []
    # assign the first Frequency partition
    fre_min_list.append(1)
    fre_max_list.append(dFreq + Frequency_Overlap)
    # assign the middle Frequency partition
    for i in range(1, Frequency_Partitions-1, 1):
        fre_min_list.append(dFreq * i - Frequency_Overlap)
        fre_max_list.append(dFreq * (i+1) + Frequency_Overlap)
    # assign the last Frequency partition
    fre_min_list.append(em_xel - (dFreq + Frequency_Overlap))
    fre_max_list.append(em_xel)

    band_params = []
    for i in range(0, Frequency_Partitions, 1):
        # convert to wavelength
        fre_min_list[i] = em_min + spectral_resolution * fre_min_list[i]
        fre_max_list[i] = em_min + spectral_resolution * fre_max_list[i]
        band_params.append(str(fre_min_list[i]) + " " + str(fre_max_list[i]))

#    return pos_params, band_params
    return range_params, band_params

def main():
    args = parseargs()
    password = casda.get_opal_password(args.opal_password, args.password_file)

    try:
        cube_vo_table = get_vo_table(args.cube_file_name, args.opal_username, password)
        # Generate dimentions of cutouts
        range_params, band_params = generate_cutouts_array(args.conf_file, cube_vo_table)

#       print(cube_vo_table['obs_publisher_did'], cube_vo_table['s_xel1'], cube_vo_table['s_xel2'], cube_vo_table['em_xel'])
#       print("POS=", pos_params)
#       print("BAND=", band_params)

        # Get access to the cube - sia call then datalink
        dataproduct_id: str = np.char.decode(cube_vo_table['obs_publisher_did'])[0]
        async_url, authenticated_id_token = casda.get_service_link_and_id(dataproduct_id,
                                                                      username = args.opal_username,
                                                                      password = password,
                                                                      destination_dir = args.destination_dir #this parameter is not optional because of a bug in casda.retrieve_direct_data_link_to_file
                                                                          )
        print (async_url, authenticated_id_token)

        # Create a job to retrieve the cutouts
        job_location = casda.create_async_soda_job([authenticated_id_token], soda_url=async_url)
        casda.add_params_to_async_job(job_location, 'POS', range_params)
        casda.add_params_to_async_job(job_location, 'BAND', band_params)
        print ('\n%d cutouts will be created' % (len(range_params)*len(band_params)))

        # Run and time the job
        run_start = time.time()
        print("\n Started running jobs: ", time.asctime(time.localtime(run_start)))
        status = casda.run_async_job(job_location)
        run_end = time.time()
        print("\n Finished: ", time.asctime(time.localtime(run_end)))
        print('\nJob finished with status %s in %.02f seconds\n\n' % (status, run_end - run_start))

        print ("Job result available at ", casda.get_results_links(job_location))

        # Optionally download
        if args.download:
           casda.download_all(job_location, args.destination_directory)
    except IOError as e:
        print(e)

    return 0

if __name__ == '__main__':
    exit(main())
