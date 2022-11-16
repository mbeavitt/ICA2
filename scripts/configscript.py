#!/bin/python3
# A simple script to write a config file. configscript.py writen on 16/11/2022 by B214618. This will not be present in the final version.

import configparser

# Initializing a configparser object
config = configparser.ConfigParser()

# Writing options to the object
config['Basic'] = {'Placeholder Option':'True'}

# Saving to file
with open('config.ini', 'w') as configfile:
    config.write(configfile)


